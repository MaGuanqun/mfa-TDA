#include "volume.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>

namespace fiber_surface {
namespace {
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
              "Boundary PLY requires IEEE 754 binary64 coordinates");

template<class T>
void order_bytes(std::array<char,sizeof(T)>& bytes) {
    const std::uint16_t one=1;
    if (*reinterpret_cast<const unsigned char*>(&one)==0)
        std::reverse(bytes.begin(),bytes.end());
}

template<class T>
void write_value(std::ostream& out, T value) {
    std::array<char,sizeof(T)> bytes;
    std::memcpy(bytes.data(),&value,sizeof(T));
    order_bytes<T>(bytes);
    out.write(bytes.data(),bytes.size());
}

template<class T>
T read_value(std::istream& in) {
    std::array<char,sizeof(T)> bytes;
    if (!in.read(bytes.data(),bytes.size())) throw std::runtime_error("Truncated boundary PLY");
    order_bytes<T>(bytes);
    T value;
    std::memcpy(&value,bytes.data(),sizeof(T));
    return value;
}
} // namespace

void write_boundary_ply(const std::string& path, const VolumeBoundary& b, bool valid) {
    std::ofstream out(path,std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write "+path);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out<<"ply\nformat binary_little_endian 1.0\n"
       <<"comment fiber_volume_boundary "<<int(valid);
    for (auto count:b.corner_components) out<<' '<<count;
    out<<' '<<b.box_contact_triangles<<"\nelement vertex "<<b.mesh.vertices.size()
       <<"\nproperty double x\nproperty double y\nproperty double z\n"
       <<"element face "<<b.mesh.triangles.size()
       <<"\nproperty list uchar int vertex_indices\nproperty uchar patch\nend_header\n";
    for (const Vec& p:b.mesh.vertices)
        for (int d=0;d<3;++d) write_value(out,p[d]);
    for (size_t i=0;i<b.mesh.triangles.size();++i) {
        write_value<std::uint8_t>(out,3);
        for (int id:b.mesh.triangles[i]) write_value<std::int32_t>(out,id);
        write_value<std::uint8_t>(out,b.patches.at(i));
    }
    out.close();
}

VolumeBoundary read_boundary_ply(const std::string& path, bool revalidate_repair) {
    std::ifstream in(path,std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open boundary PLY: "+path);
    VolumeBoundary b;
    int validated=-1;
    auto line=[&]() {
        std::string s;
        while (std::getline(in,s)) {
            if (!s.empty() && s.back()=='\r') s.pop_back();
            if (s.compare(0,8,"comment ")!=0) return s;
            std::istringstream comment(s.substr(8));
            std::string tag, extra;
            comment>>tag;
            if (tag!="fiber_volume_boundary") continue;
            if (validated!=-1 || !(comment>>validated>>b.corner_components[0]>>b.corner_components[1]
                                        >>b.corner_components[2]>>b.box_contact_triangles)
                || (validated!=0 && validated!=1) || comment>>extra)
                throw std::runtime_error("Invalid fiber boundary metadata in "+path);
        }
        throw std::runtime_error("Incomplete PLY header in "+path);
    };
    auto expect=[&](const char* declaration) {
        if (line()!=declaration) throw std::runtime_error("Expected '"+std::string(declaration)+"' in "+path);
    };
    auto count=[&](const char* name) {
        std::istringstream declaration(line());
        std::string element, type, extra;
        std::int64_t n;
        if (!(declaration>>element>>type>>n) || element!="element" || type!=name || n<0 ||
            n>std::numeric_limits<int>::max() || declaration>>extra)
            throw std::runtime_error("Invalid "+std::string(name)+" count in "+path);
        return n;
    };
    expect("ply");
    expect("format binary_little_endian 1.0");
    const auto vertices=count("vertex");
    expect("property double x");
    expect("property double y");
    expect("property double z");
    const auto faces=count("face");
    expect("property list uchar int vertex_indices");
    expect("property uchar patch");
    expect("end_header");
    if (validated<0) throw std::runtime_error("Missing fiber_volume_boundary metadata in "+path);

    // Check the fixed-width payload before allocating from the header counts.
    const auto start=in.tellg();
    in.seekg(0,std::ios::end);
    if (!in || in.tellg()-start!=vertices*24+faces*14)
        throw std::runtime_error("Boundary PLY size disagrees with its header: "+path);
    in.seekg(start);
    b.mesh.vertices.resize(vertices);
    for (Vec& p:b.mesh.vertices) {
        for (int d=0;d<3;++d) p[d]=read_value<double>(in);
        if (!p.allFinite()) throw std::runtime_error("Nonfinite vertex in "+path);
    }
    b.mesh.triangles.resize(faces);
    b.patches.reserve(faces);
    for (Triangle& triangle:b.mesh.triangles) {
        if (read_value<std::uint8_t>(in)!=3) throw std::runtime_error("Expected a triangular face in "+path);
        for (int& id:triangle) {
            id=read_value<std::int32_t>(in);
            if (id<0 || id>=vertices) throw std::runtime_error("Invalid vertex index in "+path);
        }
        const auto patch=read_value<std::uint8_t>(in);
        if (patch<1 || patch>9) throw std::runtime_error("Invalid patch ID in "+path);
        b.patches.push_back(patch);
    }
    if (!validated && !revalidate_repair)
        b.unresolved.push_back({Vec::Zero(),"The saved boundary failed validation during extraction"});
    return b;
}
} // namespace fiber_surface
