#!/bin/bash
echo "Start Runing Script"


write_vtk="./build/src/convert/write_vtk"

tracking_explicit="./build/src/critical_point_tracking/critical_point_tracking_explicit"

degenerate_case_INR="./build/src/critical_point_tracking/degenerate_case_INR"

convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"


export_raw_data="./build/src/encode/analytical/export_raw_data"


control_point_smoothing="./build/src/critical_point_tracking/control_point_smoothing"

# data_type="rotating_gaussian"
data_type="quartic_potential_2"
# data_type="sinc"
save_folder="${data_type}"






control_points="./build/src/${save_folder}/${data_type}_cpt.dat"

smoothed_control_points="./build/src/${save_folder}/${data_type}_cpt_smoothed.dat"


mfa_file="./build/src/${save_folder}/${data_type}.mfa"

smoothed_mfa_file="./build/src/${save_folder}/${data_type}_smoothed.mfa"

ori_raw_data="./build/src/${save_folder}/${data_type}_raw.dat"

step_size="16"
t_sample_ratio="16"

degenerate_point_INR="./build/src/${save_folder}/${data_type}_degenerate_INR.dat"

smoothed_degenerate_point="./build/src/${save_folder}/${data_type}_degenerate_smoothed.dat"

tracking_result="./build/src/${save_folder}/${data_type}.obj"
smoothed_tracking_result="./build/src/${save_folder}/${data_type}_smoothed.obj"


ttk_tracking_file="./build/src/${save_folder}/ttk_${data_type}.vtu"
ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt"

smoothed_ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt_smoothed.csv"


upsample_ratio="${step_size}-${step_size}-${t_sample_ratio}"

#the reshold should be really small to raw explicit function
root_finding_epsilon="1e-12"
J_threshold="1e-12"

point_itr_threshold="4.0"


if [ "${data_type}" = "quartic_potential_2" ]; then
    input_model="./CoordNet/Exp/${data_type}/super-spatial-temporal-64-5.pt"
fi



"${degenerate_case_INR}" -f "${data_type}" -b "${degenerate_point_INR}" -z "${t_sample_ratio}" -s "${step_size}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}" -m "${input_model}"


source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

# python src/python/sample_original_high_dim_func.py

# python ./src/critical_point_tracking/time_data_convert.py -i "rotating_gaussian_raw.vtk" -o "rotating_gaussian_raw.vti"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "rotating_gaussian_raw.vti" -o "rotating_gaussian_raw.csv"


# "${export_raw_data}" -d 4 -m 3 -q 4 -s 0.0 -i "${data_type}" -f "${ori_raw_data}"

# python ./src/critical_point_tracking/binary_time_data_convert.py -i "${ori_raw_data}" -o "${mfa_file}_raw.vti" -f "${data_type}"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}_raw.vti" -o "${ttk_critical_point_file}_raw.csv"

# python ./src/critical_point_tracking/time_data_convert.py -i "${mfa_file}.vtk" -o "${mfa_file}.vti"


# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}.vti" -o "${ttk_critical_point_file}.csv"



##########################################################
#smoothing and tracking on smoothed data
# gdb --args 


##########################################################
# directly work on raw explicit function

# 
# "${convert_root_to_vtk}" -f "${degenerate_point_original}" -o "${degenerate_point_original}.csv" -j 0

#gdb --args 
# "${tracking_explicit}" -f "${data_type}" -b "${tracking_result}" -z "${t_sample_ratio}" -g "${step_size}"  -x "${root_finding_epsilon}" -s "${degenerate_point_original}" -p "${point_itr_threshold}"
