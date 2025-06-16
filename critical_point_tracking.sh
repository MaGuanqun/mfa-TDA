#!/bin/bash
echo "Start Runing Script"
analytical="./build/src/encode/analytical/analytical"

gridded_3d="./build/src/encode/grid/gridded_3d"
write_vtk="./build/src/convert/write_vtk"

tracking="./build/src/critical_point_tracking/critical_point_tracking"
derivative_control_point="./build/src/critical_point/derivative_control_point"


data_type="rotating_gaussian"

save_folder="${data_type}"

control_points="./build/src/${save_folder}/${data_type}_cpt.dat"
mfa_file="./build/src/${save_folder}/${data_type}.mfa"
step_size="16"
t_sample_ratio="10"
tracking_result="./build/src/${save_folder}/${data_type}.obj"
ttk_tracking_file="./build/src/${save_folder}/ttk_${data_type}.vtu"


upsample_ratio="${step_size}-${step_size}-${t_sample_ratio}"


root_finding_epsilon="1e-10"


# if [ "${data_type}" = "rotating_gaussian" ]; then
# "${analytical}" -d 4 -m 3 -n 100 -v 28 -q 3 -s 0.0 -i "${data_type}" -f "${mfa_file}"
# else 
# "${gridded_2d}" -d 3 -f "${raw_data_file}" -i "${data_type}" -q 3 -a 0 -o "${mfa_file}"
# fi

# "${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"
# 
"${write_vtk}" -f "${mfa_file}" -t "${mfa_file}.vtk" -m 3 -d 4 -u "${upsample_ratio}" -g 0 -z 0 

# "${tracking}" -f "${mfa_file}" -b "${tracking_result}" -z "${step_size}" -a "${control_points}" -x "${root_finding_epsilon}"


source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

python ./src/critical_point_tracking/time_data_convert.py -i "${mfa_file}.vtk" -o "${mfa_file}.vti"


pvpython ./src/critical_point_tracking/tracking_script.py -i "${mfa_file}.vti" -o "${ttk_tracking_file}"