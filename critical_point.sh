#!/bin/bash
echo "Start Runing Script"
analytical="./build/src/encode/analytical/analytical"
gridded_2d="./build/src/encode/grid/gridded_2d"

derivative_control_point="./build/src/critical_point/derivative_control_point"
compute_critical_point="./build/src/critical_point/compute_critical_point"
convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"
convert_cpt_to_vtk="./build/src/critical_point/convert_cpt_to_vtk"

write_vtk="./build/src/convert/write_vtk"


data_type="sine"
data_path="./build/src/${data_type}/"
output_file_prefix="./build/src/${data_type}/"

mfa_file="${data_path}${data_type}.mfa"

raw_data_file="./build/src/ori_data/FLDSC_1_1800_3600.dat" #rti #dd07g_xxsmall_le

# "${analytical}" -d 3 -m 2 -q 3 -v 43 -n 81 -i "${data_type}" -f "${mfa_file}"

# # "${gridded_2d}" -d 3 -q 2 -f "${raw_data_file}" -i "${data_type}" -o "${mfa_file}" -a 0 -e "1.0e-2"



# "${derivative_control_point}" -f "${mfa_file}" -o "${output_file_prefix}${data_type}_cpt.dat"
# # "${convert_cpt_to_vtk}" -f "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_gradient_control"
"${compute_critical_point}" -l 1 -f "${mfa_file}" -i "${output_file_prefix}${data_type}_cpt.dat" -o "${output_file_prefix}${data_type}_cp.dat" -e 0.0999 -t 1e-10
# "${convert_root_to_vtk}" -f "${output_file_prefix}${data_type}_cp.dat" -o "${output_file_prefix}${data_type}_cp.csv" -i "${mfa_file}" -d 0


# "${write_vtk}" -f "${mfa_file}" -n "${mfa_file}.ply" -m 2 -d 3 -u 2 -g 0 -z 0 
# pvpython src/python/ttk_critical_point.py -p "${data_path}" -i "${data_type}.mfa.ply" -o "ttk.csv"