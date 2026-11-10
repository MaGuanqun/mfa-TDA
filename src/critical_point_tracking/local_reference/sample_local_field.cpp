// Sample an original-coordinate (x,y,t) box without allocating the full grid.
// MFA evaluation follows mfa_extend::recover_mfa; CoordNet follows INRModel's
// full-domain normalization and reversed (t,y,x) model input order.
#include "sampling_grid.h"
#include "../../isosurface/mfa_model.h"
#include <torch/script.h>
#include <ATen/Parallel.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <chrono>

namespace fs = std::filesystem;
using namespace local_reference;

struct Options {
    std::string backend, name, input, output, device = "cpu";
    Triple lo{}, hi{}, ratios{4,4,4}, padding{0,0,0};
    std::size_t batch_size = 65536;
    int threads = 3;
    int points_per_span = 0;
    bool region = false, dry_run = false, explicit_ratio = false;
};

void usage() {
    std::cout << "sample_local_field --backend mfa|inr --dataset NAME --input MODEL\n"
        "  --region XMIN XMAX YMIN YMAX TMIN TMAX --output PREFIX\n"
        "  [--points-per-span N | --ratio R | --ratios RX RY RT] [--padding PX PY PT]\n"
        "  [--batch-size 65536] [--threads 3] [--device cpu|cuda] [--dry-run]\n"
        "--points-per-span is MFA-only: N equal intervals per span,\n"
        "spacing = span width/N, bounds snapped outward to the domain grid.\n"
        "This mode reads MFA knots even for --dry-run. Otherwise nominal spacing\n"
        "is 0.1/ratio with exact box endpoints. Outputs: PREFIX.bin and .json.\n";
}

double number(const std::string& text) {
    std::size_t end = 0;
    const double x = std::stod(text, &end);
    if(end != text.size() || !std::isfinite(x)) throw std::invalid_argument("Invalid number: " + text);
    return x;
}

Options parse(int argc, char** argv) {
    Options o;
    auto next = [&](int& i) -> std::string {
        if(++i >= argc) throw std::invalid_argument("Missing option argument");
        return argv[i];
    };
    for(int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if(arg == "--backend") o.backend = next(i);
        else if(arg == "--dataset") o.name = next(i);
        else if(arg == "--input") o.input = next(i);
        else if(arg == "--output") o.output = next(i);
        else if(arg == "--device") o.device = next(i);
        else if(arg == "--dry-run") o.dry_run = true;
        else if(arg == "--region") {
            for(int a = 0; a < 3; ++a) {o.lo[a] = number(next(i)); o.hi[a] = number(next(i));}
            o.region = true;
        } else if(arg == "--ratio") {o.ratios.fill(number(next(i))); o.explicit_ratio = true;}
        else if(arg == "--ratios" || arg == "--padding") {
            auto& target = arg == "--ratios" ? o.ratios : o.padding;
            for(auto& x : target) x = number(next(i));
            if(arg == "--ratios") o.explicit_ratio = true;
        } else if(arg == "--batch-size" || arg == "--threads" || arg == "--points-per-span") {
            double x = number(next(i));
            if(x < 1 || x > std::numeric_limits<int>::max() || std::floor(x) != x)
                throw std::invalid_argument(arg + " must be a positive integer <= INT_MAX");
            if(arg == "--points-per-span") {
                if(x < 2) throw std::invalid_argument("--points-per-span must be at least 2 (TTK needs three frames)");
                o.points_per_span = static_cast<int>(x);
            } else if(arg == "--threads") o.threads = static_cast<int>(x);
            else o.batch_size = static_cast<std::size_t>(x);
        } else throw std::invalid_argument("Unknown option: " + arg);
    }
    if(o.backend != "mfa" && o.backend != "inr") throw std::invalid_argument("--backend must be mfa or inr");
    if(!o.region || o.output.empty() || o.input.empty())
        throw std::invalid_argument("--region, --input and --output are required");
    if(o.device != "cpu" && o.device != "cuda") throw std::invalid_argument("--device must be cpu or cuda");
    if(o.points_per_span && (o.backend != "mfa" || o.explicit_ratio))
        throw std::invalid_argument("--points-per-span requires --backend mfa and cannot be combined with --ratio/--ratios");
    return o;
}

std::string json_string(const std::string& s) {
    std::ostringstream out; out << '"';
    for(unsigned char c : s) {
        if(c == '"' || c == '\\') out << '\\' << c;
        else if(c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    return out.str() + '"';
}

template<class T> void json_array(std::ostream& out, const std::array<T,3>& a) {
    out << '[' << a[0] << ',' << a[1] << ',' << a[2] << ']';
}

void metadata(std::ostream& out, const Options& o, const Dataset& d, const Grid& g) {
    out << std::setprecision(17) << "{\n  \"format\": \"local-field-v1\",\n"
        << "  \"backend\": " << json_string(o.backend) << ",\n  \"dataset\": " << json_string(o.name)
        << ",\n  \"model\": " << json_string(fs::absolute(o.input).string())
        << ",\n  \"dtype\": \"<f8\",\n  \"layout\": \"t,y,x\",\n  \"dimensions\": "; json_array(out,g.dimensions);
    out << ",\n  \"domain_min\": "; json_array(out,d.min);
    out << ",\n  \"domain_max\": "; json_array(out,d.max);
    out << ",\n  \"region_min\": "; json_array(out,o.lo);
    out << ",\n  \"region_max\": "; json_array(out,o.hi);
    out << ",\n  \"sample_min\": "; json_array(out,g.min);
    out << ",\n  \"sample_max\": "; json_array(out,g.max);
    out << ",\n  \"spacing\": "; json_array(out,g.spacing);
    out << ",\n  \"sampling_mode\": " << json_string(o.points_per_span ? "mfa-knot-spans" : "physical-ratio");
    out << ",\n  \"points_per_span\": ";
    if(o.points_per_span) out << o.points_per_span; else out << "null";
    out << ",\n  \"mfa_span_counts\": ";
    if(o.points_per_span) json_array(out,d.base_intervals); else out << "null";
    out << ",\n  \"ratios\": ";
    if(o.points_per_span) out << "null"; else json_array(out,o.ratios);
    out << ",\n  \"padding\": "; json_array(out,o.padding);
    out << "\n}\n";
}

int main(int argc, char** argv) {
    if(argc == 1 || (argc == 2 && std::string(argv[1]) == "--help")) {usage(); return 0;}
    try {
        auto o = parse(argc, argv);
        auto d = dataset(o.name);
        const unsigned short endian = 1;
        if(*reinterpret_cast<const unsigned char*>(&endian) != 1)
            throw std::runtime_error("This binary writer requires a little-endian host");
        if((!o.dry_run || o.points_per_span) && !fs::is_regular_file(o.input))
            throw std::runtime_error("Model does not exist: " + o.input);
        tbb::global_control limit(tbb::global_control::max_allowed_parallelism, o.threads);
        at::set_num_threads(o.threads);
        torch::NoGradGuard no_grad;
        // MPI is only needed for DIY model I/O, not for CoordNet sampling.
        std::unique_ptr<diy::mpi::environment> mpi;
        std::shared_ptr<marching_triangles::MfaModel> mfa_model;
        torch::jit::Module module;
        const torch::Device device(o.device);
        auto dtype = torch::kFloat64;
        if(o.backend == "mfa" && (!o.dry_run || o.points_per_span)) {
            mpi = std::make_unique<diy::mpi::environment>(argc,argv);
            diy::mpi::communicator world;
            mfa_model = std::make_shared<marching_triangles::MfaModel>(o.input,world);
            const auto& b = mfa_model->block();
            if(!b.mfa || b.mfa->dom_dim != 3 || b.mfa->nvars() != 1 || b.mfa->var(0).dim() != 1)
                throw std::runtime_error("Expected one scalar MFA field with domain (x,y,t)");
            if(o.points_per_span) {
                const auto& mesh = b.mfa->var(0).tmesh;
                if(mesh.tensor_prods.size() != 1 || mesh.all_knots.size() != 3)
                    throw std::runtime_error("--points-per-span requires a single tensor-product MFA field");
                for(int a = 0; a < 3; ++a) {
                    if(std::abs(b.core_mins[a]-d.min[a]) > 1e-10 || std::abs(b.core_maxs[a]-d.max[a]) > 1e-10)
                        throw std::runtime_error("MFA core domain does not match the selected dataset");
                    d.base_intervals[a] = uniform_span_count(mesh.all_knots[a]);
                }
                o.ratios.fill(o.points_per_span);
            }
        } else if(o.backend == "inr" && !o.dry_run) {
            module = torch::jit::load(o.input, device);
            module.eval();
            for(const auto& p : module.parameters()) {dtype = p.scalar_type(); break;}
            if(dtype != torch::kFloat32 && dtype != torch::kFloat64)
                throw std::runtime_error("CoordNet parameters must be float32 or float64");
        }
        const Grid g(d, o.lo, o.hi, o.ratios, o.padding, o.points_per_span != 0);
        if(mfa_model) {
            const auto& b = mfa_model->block();
            for(int a = 0; a < 3; ++a)
                if(!std::isfinite(b.core_mins[a]) || !std::isfinite(b.core_maxs[a]) ||
                   b.core_mins[a] >= b.core_maxs[a] || g.min[a] < b.core_mins[a] || g.max[a] > b.core_maxs[a])
                    throw std::runtime_error("Requested sample region lies outside the MFA model domain");
        }
        std::cout << "Grid (nx,ny,nt): " << g.dimensions[0] << ',' << g.dimensions[1] << ',' << g.dimensions[2]
                  << "; scalar storage: " << g.size * 8. / (1024.*1024.*1024.) << " GiB\n";
        if(o.points_per_span) {
            std::cout << "MFA spans: "; json_array(std::cout,d.base_intervals);
            std::cout << "; spacing: "; json_array(std::cout,g.spacing); std::cout << '\n';
        }
        if(o.dry_run) {metadata(std::cout,o,d,g); return 0;}
        const fs::path prefix(o.output);
        if(!prefix.parent_path().empty()) fs::create_directories(prefix.parent_path());
        const std::string temporary = o.output + ".bin.tmp";
        std::ofstream output(temporary,std::ios::binary);
        output.exceptions(std::ios::failbit | std::ios::badbit);
        auto last_progress = std::chrono::steady_clock::now();
        for(std::size_t start = 0; start < g.size;) {
            const auto count = std::min(o.batch_size,g.size-start);
            std::vector<double> values(count);
            if(mfa_model) {
                const auto& b = mfa_model->block();
                // DecodeVar constructs a decoder and scratch buffers per point.
                // Use its same VolPt evaluation with buffers reused per task.
                tbb::parallel_for(tbb::blocked_range<std::size_t>(0,count,1024),[&](const auto& range) {
                    const auto& field = b.mfa->var(0);
                    mfa::Decoder<double> decoder(field,0);
                    mfa::DecodeInfo<double> info(field,Eigen::VectorXi());
                    Eigen::VectorXd param(3), value(1);
                    for(auto i = range.begin(); i != range.end(); ++i) {
                        const auto p = g.point(start+i);
                        for(int a = 0; a < 3; ++a)
                            param[a] = (p[a]-b.core_mins[a])/(b.core_maxs[a]-b.core_mins[a]);
                        decoder.VolPt(param,value,info,field.tmesh.tensor_prods[0]);
                        values[i] = value[0];
                    }
                });
            } else {
                auto input = torch::empty({static_cast<int64_t>(count),3},torch::kFloat64);
                auto coords = input.accessor<double,2>();
                for(std::size_t i = 0; i < count; ++i) {
                    const auto p = g.point(start+i);
                    for(int a = 0; a < 3; ++a)
                        coords[i][2-a] = 2*(p[a]-d.min[a])/(d.max[a]-d.min[a])-1;
                }
                auto result = module.forward({input.to(device,dtype)}).toTensor()
                    .to(torch::kCPU,torch::kFloat64).contiguous();
                if(static_cast<std::size_t>(result.numel()) != count)
                    throw std::runtime_error("CoordNet must output one scalar per input point");
                std::copy_n(result.data_ptr<double>(),count,values.data());
            }
            for(double value : values) if(!std::isfinite(value)) throw std::runtime_error("Nonfinite sampled field value");
            output.write(reinterpret_cast<const char*>(values.data()),count*sizeof(double));
            start += count;
            const auto now = std::chrono::steady_clock::now();
            if(now-last_progress > std::chrono::seconds(30) || start == g.size) {
                std::cout << "Sampled " << start << '/' << g.size << " points ("
                          << 100.*start/g.size << "%)" << std::endl;
                last_progress = now;
            }
        }
        output.close();
        std::ofstream info(o.output+".json.tmp");
        info.exceptions(std::ios::failbit | std::ios::badbit);
        metadata(info,o,d,g); info.close();
        fs::rename(temporary,o.output+".bin");
        fs::rename(o.output+".json.tmp",o.output+".json");
        std::cout << "Wrote " << o.output << ".bin and .json\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "sample_local_field: " << e.what() << '\n'; return 1;
    }
}
