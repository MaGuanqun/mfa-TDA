#!/bin/bash
echo "Start Runing Script"
analytical="./build/src/encode/analytical/analytical"
gridded_2d="./build/src/encode/grid/gridded_2d"

jacobi_set="./build/src/contour/jacobi_set"
write_two_value="./build/src/convert/write_two_value"

convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"
write_vtk="./build/src/convert/write_vtk"

derivative_control_point="./build/src/critical_point/derivative_control_point"
compute_critical_point="./build/src/critical_point/compute_critical_point"
compute_critical_point_jacobi_set="./build/src/contour/compute_critical_point_jacobi_set"


extract_ttk_js="./src/python/extract_jacobi_set.py"

count_betti_num="./src/contour/count_betti_num.py"

# data_type="hurricane_isabel"
data_type="vortex_street"
# data_type="gaussian_mixture"
# data_type="boussinesq"
save_folder="${data_type}"

extension="ply"

mfa_file="./build/src/${save_folder}/${data_type}.mfa"
mfa_file2="./build/src/${save_folder}/${data_type}_2.mfa"
root_file="./build/src/${save_folder}/${data_type}_root.dat"
root_file_csv="./build/src/${save_folder}/${data_type}_root.csv"


ridge_valley_file1="./build/src/${save_folder}/${data_type}_ridge_valley_graph"
ridge_valley_file_csv1="./build/src/${save_folder}/${data_type}_ridge_valley_graph.csv"


ridge_valley_file2="./build/src/${save_folder}/${data_type}_ridge_valley_graph2.obj"
ridge_valley_file_csv2="./build/src/${save_folder}/${data_type}_ridge_valley_graph2.csv"


control_points="./build/src/${save_folder}/${data_type}_cpt.dat"
control_points2="./build/src/${save_folder}/${data_type}_cpt2.dat"

critical_point_jacobi_set="./build/src/${save_folder}/${data_type}_rt_jacobi_set.dat"
critical_point_jacobi_set_csv_file="./build/src/${save_folder}/${data_type}_rt_jacobi_set.csv"

critical_point_file="./build/src/${save_folder}/${data_type}_rt.dat"
critical_point_file2="./build/src/${save_folder}/${data_type}_rt2.dat"
critical_point_csv_file="./build/src/${save_folder}/${data_type}_rt.csv"
critical_point_csv_file2="./build/src/${save_folder}/${data_type}_rt2.csv"


step_size="2"
step_size2="32"
step_size3="8"
step_size4="16"
step_size5="32"
step_size6="64"
step_size7="128"
step_size8="256"



root_finding_epsilon2="1e-6"
root_finding_epsilon3="1e-8"
root_finding_epsilon4="1e-12"


connection_threshold="2.0"

connection_threshold2="1.0"
connection_threshold3="1.5"
connection_threshold4="2.5"



if [ "${data_type}" = "boussinesq" ]; then
    raw_data_file="./build/src/boussinesq/velocity_t2000.bin"  #    
    raw_data_file2="./build/src/boussinesq/velocity_t2001.bin" 
    root_finding_epsilon="1e-10"
    trace_split_grad_square_threshold="1e-4"
elif [ "${data_type}" = "hurricane_isabel" ]; then
    raw_data_file="./build/src/hurricane_isabel/Pf30_50.bin"  #    
    raw_data_file2="./build/src/hurricane_isabel/TCf30_50.bin" 
    root_finding_epsilon="1e-10"
    trace_split_grad_square_threshold="1e-4"
elif [ "${data_type}" = "vortex_street" ]; then
    raw_data_file="./build/src/vortex_street/velocity_t1500.bin" 
    raw_data_file2="./build/src/vortex_street/velocity_t1501.bin" 
    root_finding_epsilon="1e-10"
    trace_split_grad_square_threshold="1e-4"
elif [ "${data_type}" = "gaussian_mixture" ]; then
    root_finding_epsilon="1e-10"
    trace_split_grad_square_threshold="1e-4"
fi



if [ "${data_type}" = "gaussian_mixture" ]; then
    "${analytical}" -d 3 -m 2 -n 101 -v 21 -q 4 -s 0.0 -i "${data_type}1" -f "${mfa_file}"
    "${analytical}" -d 3 -m 2 -n 101 -v 21 -q 4 -s 0.0 -i "${data_type}2" -f "${mfa_file2}"
else
    "${gridded_2d}" -d 3 -f "${raw_data_file}" -i "${data_type}" -q 4 -a 0 -o "${mfa_file}"
    "${gridded_2d}" -d 3 -f "${raw_data_file2}" -i "${data_type}" -q 4 -a 0 -o "${mfa_file2}"
fi


"${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"
"${derivative_control_point}" -f "${mfa_file2}" -o "${control_points2}"

same_root_epsilon=$(awk "BEGIN {print 1 / ${step_size2}}")
echo "Same root epsilon: ${same_root_epsilon}"

"${compute_critical_point_jacobi_set}" -l 1 -f "${mfa_file}" -g "${mfa_file2}" -i "${control_points}" -j "${control_points2}" -o "${critical_point_jacobi_set}" -v "1e-6" -t 1e-10 -e "${same_root_epsilon}"


# "${convert_root_to_vtk}" -i "${mfa_file}" -o "${critical_point_jacobi_set_csv_file}" -f "${critical_point_jacobi_set}" -e "${same_root_epsilon}"



# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size}.obj" -z "${step_size}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

"${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -b "${ridge_valley_file1}${step_size2}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon2}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon3}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon4}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1


# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold2}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1 -a "${connection_threshold2}"

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold3}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1 -a "${connection_threshold3}"

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size2}_${connection_threshold4}.obj" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1 -a "${connection_threshold4}"



# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size3}.obj" -z "${step_size3}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size4}.obj" -z "${step_size4}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size5}.obj" -z "${step_size5}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size6}.obj" -z "${step_size6}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size7}.obj" -z "${step_size7}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file1}${step_size8}.obj" -z "${step_size8}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# "${jacobi_set}" -f "${mfa_file}" -g "${mfa_file2}" -r "${root_file}" -b "${ridge_valley_file2}" -z "${step_size2}" -v 0 -p "${trace_split_grad_square_threshold}" -y 0.5 -x "${root_finding_epsilon}" -j "${critical_point_jacobi_set}" -c "${critical_point_jacobi_set}" -o 1

# 
# gdb "${jacobi_set}" -ex "set args -f ${mfa_file} -g ${mfa_file2} -r ${root_file} -b ${ridge_valley_file1} -z ${step_size} -v 0 -p ${trace_split_grad_square_threshold} -y 0.5 -x ${root_finding_epsilon} -j ${critical_point_jacobi_set} -c ${critical_point_jacobi_set}"




# "${write_two_value}" -f "${mfa_file}" -g "${mfa_file2}" -t "${mfa_file}.vtk" -m 2 -d 3 -u "${step_size}"
"${write_two_value}" -f "${mfa_file}" -g "${mfa_file2}" -t "${mfa_file}_up.vtk" -m 2 -d 3 -u "${step_size2}"


# 4.5589 for ackley


# "${analytical}" -d 3 -m 2 -n 401 -v 101 -q 4 -s 0.0 -i "gaussian_pair" -o "${mfa_file}"

# "${write_h}" -f "${mfa_file}" -t "${mfa_file}_h.vtk" -m 2 -d 3 -u 10 -s "0.85-1.0-0.85-1.0"


source ~/enter/etc/profile.d/conda.sh
conda activate mfa

# python "${extract_ttk_js}" "${mfa_file}.vtk" "./build/src/${save_folder}/js_ttk_${step_size}.obj" "./build/src/${save_folder}/js_ttk_${step_size}.bin"

python "${extract_ttk_js}" "${mfa_file}_up.vtk" "./build/src/${save_folder}/js_ttk_${step_size2}.obj" "./build/src/${save_folder}/js_ttk_${step_size2}.bin"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size}.obj"
python "${count_betti_num}" "${ridge_valley_file1}${step_size2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${root_finding_epsilon4}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${connection_threshold2}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${connection_threshold3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size2}_${connection_threshold4}.obj"

# python "${count_betti_num}" "${ridge_valley_file1}${step_size3}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size4}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size5}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size6}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size7}.obj"
# python "${count_betti_num}" "${ridge_valley_file1}${step_size8}.obj"
# python "${count_betti_num}" "${ridge_valley_file2}"

# python "${count_betti_num}" "./build/src/${save_folder}/js_ttk_${step_size}.obj"
python "${count_betti_num}" "./build/src/${save_folder}/js_ttk_${step_size2}.obj"
# 
# ttk_contour_error_file="./build/src/contour/compute_function_error"
# "${ttk_contour_error_file}" -f "${mfa_file}" -g "${mfa_file2}"  -r "./build/src/${save_folder}/js_ttk_${step_size}.bin" -s "./build/src/${save_folder}/js_ttk_${step_size}_error.dat" -z "2"
# # "${ttk_contour_error_file}" -f "${mfa_file}" -g "${mfa_file2}"  -r "./build/src/${save_folder}/js_ttk_${step_size2}.bin" -s "./build/src/${save_folder}/js_ttk_${step_size}_error2.dat" -z "2"