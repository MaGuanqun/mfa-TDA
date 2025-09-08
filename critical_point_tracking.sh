#!/bin/bash
echo "Start Runing Script"
analytical="./build/src/encode/analytical/analytical"

gridded_3d="./build/src/encode/grid/gridded_3d"
write_vtk="./build/src/convert/write_vtk"

tracking="./build/src/critical_point_tracking/critical_point_tracking"
derivative_control_point="./build/src/critical_point/derivative_control_point"

degenerate_case="./build/src/critical_point_tracking/degenerate_case"
convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"


export_raw_data="./build/src/encode/analytical/export_raw_data"

control_point_smoothing="./build/src/critical_point_tracking/control_point_smoothing"

data_type="rotating_gaussian"
# data_type="sinc"
save_folder="${data_type}"






control_points="./build/src/${save_folder}/${data_type}_cpt.dat"

smoothed_control_points="./build/src/${save_folder}/${data_type}_cpt_smoothed.dat"


mfa_file="./build/src/${save_folder}/${data_type}.mfa"

smoothed_mfa_file="./build/src/${save_folder}/${data_type}_smoothed.mfa"

ori_raw_data="./build/src/${save_folder}/${data_type}_raw.dat"

step_size="16"
t_sample_ratio="16"

degenerate_point="./build/src/${save_folder}/${data_type}_degenerate.dat"

smoothed_degenerate_point="./build/src/${save_folder}/${data_type}_degenerate_smoothed.dat"

tracking_result="./build/src/${save_folder}/${data_type}.obj"
smoothed_tracking_result="./build/src/${save_folder}/${data_type}_smoothed.obj"


ttk_tracking_file="./build/src/${save_folder}/ttk_${data_type}.vtu"
ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt"

smoothed_ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt_smoothed.csv"


upsample_ratio="${step_size}-${step_size}-${t_sample_ratio}"


root_finding_epsilon="1e-10"
J_threshold="1e-10"

point_itr_threshold="4.0"

# block="0.45-0.55-0.6225-0.7225-0.21863-0.31863"

# block="0.4449-0.5449-0.2762-0.3762-0.6809-0.7809"

# if [ "${data_type}" = "rotating_gaussian" ]; then
# "${analytical}" -d 4 -m 3 -q 3 -s 0.0 -i "${data_type}" -f "${mfa_file}"
# else 
# "${gridded_2d}" -d 3 -f "${raw_data_file}" -i "${data_type}" -q 3 -a 0 -o "${mfa_file}"
# fi

# "${export_raw_data}" -d 4 -m 3 -q 4 -s 0.0 -i "${data_type}" -f "${ori_raw_data}"


# "${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"



# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}.vtk" -m 3 -d 4 -u "${upsample_ratio}" -g 0 -z 0 #-s "${block}"

"${degenerate_case}" -f "${mfa_file}" -b "${degenerate_point}" -z "${t_sample_ratio}" -s "${step_size}" -a "${control_points}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}"

"${convert_root_to_vtk}" -f "${degenerate_point}" -o "${degenerate_point}.csv" -i "${mfa_file}" -d 0



# "${tracking}" -f "${mfa_file}" -b "${tracking_result}" -z "${t_sample_ratio}" -g "${step_size}"  -a "${control_points}" -x "${root_finding_epsilon}" -s "${degenerate_point}" -p "${point_itr_threshold}"




source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

# python src/python/sample_original_high_dim_func.py

# python ./src/critical_point_tracking/time_data_convert.py -i "rotating_gaussian_raw.vtk" -o "rotating_gaussian_raw.vti"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "rotating_gaussian_raw.vti" -o "rotating_gaussian_raw.csv"

# python ./src/critical_point_tracking/binary_time_data_convert.py -i "${ori_raw_data}" -o "${mfa_file}_raw.vti" -f "${data_type}"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}_raw.vti" -o "${ttk_critical_point_file}_raw.csv"

# python ./src/critical_point_tracking/time_data_convert.py -i "${mfa_file}.vtk" -o "${mfa_file}.vti"


# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}.vti" -o "${ttk_critical_point_file}.csv"


# gdb --args 
# "${control_point_smoothing}" -f "${mfa_file}" -o "${smoothed_mfa_file}" -s 0.5

# "${derivative_control_point}" -f "${smoothed_mfa_file}" -o "${smoothed_control_points}"

# "${degenerate_case}" -f "${smoothed_mfa_file}" -b "${smoothed_degenerate_point}" -z "${t_sample_ratio}" -s "${step_size}" -a "${smoothed_control_points}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}"

# "${convert_root_to_vtk}" -f "${smoothed_degenerate_point}" -o "${smoothed_degenerate_point}.csv" -i "${smoothed_mfa_file}" -d 0

# "${tracking}" -f "${smoothed_mfa_file}" -b "${smoothed_tracking_result}" -z "${t_sample_ratio}" -g "${step_size}"  -a "${smoothed_control_points}" -x "${root_finding_epsilon}" -s "${smoothed_degenerate_point}" -p "${point_itr_threshold}"


# "${write_vtk}" -f "${smoothed_mfa_file}" -t "${smoothed_mfa_file}.vtk" -g 0 -z 0 -i "${data_type}" -u "${upsample_ratio}" #-s "${block}"

# python ./src/critical_point_tracking/time_data_convert.py -i "${smoothed_mfa_file}.vtk" -o "${smoothed_mfa_file}.vti"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${smoothed_mfa_file}.vti" -o "${smoothed_ttk_critical_point_file}.csv"

# pvpython ./src/critical_point_tracking/tracking_script.py -i "${mfa_file}.vti" -o "${ttk_tracking_file}"