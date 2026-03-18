#!/bin/bash
echo "Start Runing Script"
analytical="./build/src/encode/analytical/analytical"

time_varying="./build/src/encode/timeseries/time_varying"
write_vtk="./build/src/convert/write_vtk"
write_slice_bin="./build/src/convert/write_bin_select_slices"

tracking="./build/src/critical_point_tracking/critical_point_tracking"
tracking_explicit="./build/src/critical_point_tracking/critical_point_tracking_explicit"

count_betti_num="./src/contour/count_betti_num.py"

derivative_control_point="./build/src/critical_point/derivative_control_point"

degenerate_case="./build/src/critical_point_tracking/degenerate_case"
degenerate_case_explicit="./build/src/critical_point_tracking/degenerate_case_explicit"

convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"


export_raw_data="./build/src/encode/analytical/export_raw_data"


control_point_smoothing="./build/src/critical_point_tracking/control_point_smoothing"

data_type="vortex" #quartic_potential_3d vortex

# data_type="quartic_potential_2"

# data_type="sinc"
save_folder="${data_type}"



convert_obj_domain=./src/critical_point_tracking/convert_obj_back_to_domain.py


control_points="./build/src/${save_folder}/${data_type}_cpt.dat"

smoothed_control_points="./build/src/${save_folder}/${data_type}_cpt_smoothed.dat"


mfa_file="./build/src/${save_folder}/${data_type}.mfa"

smoothed_mfa_file="./build/src/${save_folder}/${data_type}_smoothed.mfa"

ori_raw_data="./build/src/${save_folder}/${data_type}_raw.dat"



degenerate_point="./build/src/${save_folder}/${data_type}_degenerate"
degenerate_point_original="./build/src/${save_folder}/${data_type}_degenerate"

smoothed_degenerate_point="./build/src/${save_folder}/${data_type}_degenerate_smoothed.dat"

tracking_result="./build/src/${save_folder}/${data_type}"
tracking_result_w_type="./build/src/${save_folder}/${data_type}"


INR_tracking_result="./build/src/${save_folder}_INR/${data_type}"


smoothed_tracking_result="./build/src/${save_folder}/${data_type}_smoothed.obj"

ori_tracking_result="./build/src/${save_folder}/ori_${data_type}"
edge_type_file="./build/src/${save_folder}/${data_type}_edge_type"
ori_tracking_result_w_type="./build/src/${save_folder}/${data_type}_original"

ttk_tracking_file="./build/src/${save_folder}/ttk_${data_type}.vtu"
ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt"

smoothed_ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt_smoothed.csv"


boundary_point_INR="./build/src/${save_folder}/${data_type}_boundary"

filter_points="./src/python/filter_point.py"



#the reshold should be really small to raw explicit function
root_finding_epsilon="1e-10"
J_threshold="1e-6"

point_itr_threshold="4.0"




raw_data_file="./build/src/ori_data/vorts/vorts_list.txt"



# block="0.45-0.55-0.6225-0.7225-0.21863-0.31863"

# block="0.4449-0.5449-0.2762-0.3762-0.6809-0.7809"

# if [ "${data_type}" = "rotating_gaussian" ]; then
# "${analytical}" -d 4 -m 3 -q 3 -s 0.0 -i "${data_type}" -f "${mfa_file}"
# else 
# "${time_varying}" -f "${raw_data_file}" -i "${data_type}" -q 3 -o "${mfa_file}" -s 0
# fi


# "${derivative_control_point}" -f "${mfa_file}" -o "${control_points}"


# for step_size in 3 2
# do
# # t_sample_ratio="64"


# step_size=96
# for J_threshold in "1e-2" # "1e-3" "1e-4" "1e-5" "1e-6" "1e-7" "1e-8"
# do  
# echo "degenerate case with threshold: ${J_threshold}"
# "${degenerate_case}" -f "${mfa_file}" -b "${degenerate_point}" -z "${step_size}" -s "${step_size}" -a "${control_points}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}"
# done
# done


for step_size in 8 16
do
    # "${convert_root_to_vtk}" -f "${degenerate_point}${step_size}.dat" -o "${degenerate_point}_${step_size}.csv" -j 0 -t 0
    # # gdb --args 
    # for root_finding_epsilon in "1e-6" #"1e-7" "1e-8" "1e-9" "1e-10" "1e-11" "1e-12" "1e-13"
    # do
    echo "test"
    # echo "tracking with step size: ${step_size}"
    "${tracking}" -f "${mfa_file}" -b "${tracking_result}_${step_size}.obj" -z "${step_size}" -g "${step_size}"  -a "${control_points}" -x "${root_finding_epsilon}" -s "${degenerate_point}${step_size}.dat" -p "${point_itr_threshold}" -e "${edge_type_file}_${step_size}.csv" -j "${boundary_point_INR}"
    
    
    # "${convert_root_to_vtk}" -f "${boundary_point_INR}${step_size}.dat" -o "${boundary_point_INR}${step_size}.csv" -j 0 -t 0

    # python ./src/python/remove_obstacle.py --data "${data_type}" --input_csv "${boundary_point_INR}${step_size}.csv" --output_csv "${boundary_point_INR}${step_size}_remove_obstacle.csv" --for_our_csv 1


    # python "${filter_points}" -i "${boundary_point_INR}${step_size}_remove_obstacle.csv"
    # source ~/enter/etc/profile.d/conda.sh
    # conda activate mfa
    # python $convert_obj_domain --csv "${degenerate_point}_${step_size}.csv" --function "${data_type}" --obj "${tracking_result}_${step_size}.obj"
    # python ./src/python/merge_obj_edge_type.py -i "${tracking_result}_${step_size}.obj" -j "${edge_type_file}_${step_size}.csv" -o "${tracking_result}_${step_size}.vtp"

    # pvpython "${count_betti_num}" "${tracking_result}_${step_size}.vtp"

    # upsample_ratio="${step_size}-${step_size}-${step_size}"


    # pvpython ./src/python/remove_obstacle.py -i "${tracking_result}_${step_size}.vtp" -o "${tracking_result}_${step_size}_remove_obstacle.vtp" --data "${data_type}"
    # #     echo "remove obstacle/n"
    # pvpython "${count_betti_num}" "${tracking_result}_${step_size}_remove_obstacle.vtp"

    # "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}_${step_size}.bin" -m 3 -d 4 -u "${upsample_ratio}" -g 0 -z 0 -b 1 #-s "${block}"
    # python ./src/critical_point_tracking/binary_time_data_convert_mfa.py -i "${mfa_file}_${step_size}.bin" -o "${mfa_file}_${step_size}.vti" -f "${data_type}" --float_type float64 --step_size "${step_size}"
    # source ~/enter/etc/profile.d/conda.sh
    # conda activate mfa_env
    # python ./src/critical_point_tracking/time_data_convert.py -i "${mfa_file}_${step_size}.vtk" -o "${mfa_file}_${step_size}.vti"
    # pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}_${step_size}.vti" -o "${mfa_file}_${step_size}.csv"

    # pvpython ./src/critical_point_tracking/mfa_bin_critical_point.py --cpp_exe "${write_slice_bin}" --chunk_size 20 --step_size "${step_size}" --float_type 'float64' --function "${data_type}" --server 0 --output_csv "${mfa_file}_${step_size}.csv" --mfa_file "${mfa_file}"

    done
# done

# for step_size in "2" "4" "8" "16" "32" "64"
# do

# done



# for step_size in 16
# do
# #step size
#     if [ "${step_size}" = "16" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.00625" 
#             t_step_size="0.00520833"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00694444"
#             t_step_size="0.0133929"
#         fi
#     elif [ "${step_size}" = "8" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.0125" 
#             t_step_size="0.01041666"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.01388888"
#             t_step_size="0.0267858"
#         fi
#     elif [ "${step_size}" = "2" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.05" 
#             t_step_size="0.04166666"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.05555555"
#             t_step_size="0.1071432"
#         fi
#     elif [ "${step_size}" = "3" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.03333333" 
#             t_step_size="0.02777777"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.03703703"
#             t_step_size="0.0714283"
#         fi
#     elif [ "${step_size}" = "6" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.01666666" 
#             t_step_size="0.01388888"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.01851851"
#             t_step_size="0.03571415"
#         fi
#     elif [ "${step_size}" = "12" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.00833333" 
#             t_step_size="0.00694444"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00925925"
#             t_step_size="0.017857075"
#         fi
#     elif [ "${step_size}" = "24" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.00416666" 
#             t_step_size="0.00347222"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00462962"
#             t_step_size="0.0089285375"
#         fi
#     elif [ "${step_size}" = "48" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.00208333" 
#             t_step_size="0.00173611"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00231481"
#             t_step_size="0.00446426875"
#         fi
#     elif [ "${step_size}" = "96" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.001041666" 
#             t_step_size="0.0008680555"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.001157407"
#             t_step_size="0.002232134375"
#         fi
#     elif [ "${step_size}" = "4" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.025" 
#             t_step_size="0.02083333"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.02777777"
#             t_step_size="0.0535716"
#         fi
#     elif [ "${step_size}" = "32" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.003125" 
#             t_step_size="0.002604166"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00347222"
#             t_step_size="0.00669645"
#         fi
#     elif [ "${step_size}" = "64" ]; then
#         if [ "${data_type}" = "vortex_street_3d" ]; then
#             spatial_step_size="0.0015625" 
#             t_step_size="0.0013020833"
#         elif [ "${data_type}" = "boussinesq_3d" ]; then
#             spatial_step_size="0.00173611"
#             t_step_size="0.003348225"
#         fi
#     fi
#     # echo "step size $step_size"

# # source ~/enter/etc/profile.d/conda.sh
# # conda activate siren
# # python ./src/critical_point_tracking/matching_ratio.py "${tracking_result}_${step_size}.vtp" "${mfa_file}_${step_size}.csv" "${spatial_step_size}" "${t_step_size}" --header

# # python ./src/python/remove_obstacle.py --data "${data_type}" --input_csv "${mfa_file}_${step_size}.csv" --output_csv "${mfa_file}_${step_size}_remove_obstacle.csv"

# # python ./src/critical_point_tracking/matching_ratio.py "${tracking_result}_${step_size}_remove_obstacle.vtp" "${mfa_file}_${step_size}_remove_obstacle.csv" "${spatial_step_size}" "${t_step_size}" --header

# # INR_tracking_result="./build/src/${save_folder}_INR/${data_type}"

# # spatial_step_size2=0.0125
# # t_step_size2=0.0125

# # python ./src/critical_point_tracking/matching_ratio.py "${tracking_result}_${step_size}.vtp" "${INR_tracking_result}_8.vtp" "${spatial_step_size2}" "${t_step_size2}" --header


# done
# if [ "${data_type}" = "vortex_street_3d" ]; then
#     spatial_step_size="0.0125" 
#     t_step_size="0.01041666"
# elif [ "${data_type}" = "boussinesq_3d" ]; then
#     spatial_step_size="0.01388888"
#     t_step_size="0.0220588"
# fi





# source ~/enter/etc/profile.d/conda.sh
# conda activate mfa_env

# python src/python/sample_original_high_dim_func.py --function_name "${data_type}" --output_name "${ori_raw_data}.vtk"

# python ./src/critical_point_tracking/time_data_convert.py -i "${ori_raw_data}.vtk" -o "${ori_raw_data}.vti"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "rotating_gaussian_raw.vti" -o "rotating_gaussian_raw.csv"


# "${export_raw_data}" -d 4 -m 3 -q 4 -s 0.0 -i "${data_type}" -f "${ori_raw_data}"

# python ./src/critical_point_tracking/binary_time_data_convert.py -i "${ori_raw_data}" -o "${mfa_file}_raw.vti" -f "${data_type}"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}_raw.vti" -o "${ttk_critical_point_file}_raw.csv"

# "${write_vtk}" -f "${mfa_file}" -t "${mfa_file}.vtk" -m 3 -d 4 -u "${upsample_ratio}" -g 0 -z 0 #-s "${block}"
# python ./src/critical_point_tracking/time_data_convert.py -i "${mfa_file}.vtk" -o "${mfa_file}.vti"

# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${mfa_file}.vti" -o "${ttk_critical_point_file}_${step_size}.csv"



##########################################################
#smoothing and tracking on smoothed data
# gdb --args 
# "${control_point_smoothing}" -f "${mfa_file}" -o "${smoothed_mfa_file}" -s 0.5

# "${derivative_control_point}" -f "${smoothed_mfa_file}" -o "${smoothed_control_points}"

# "${degenerate_case}" -f "${smoothed_mfa_file}" -b "${smoothed_degenerate_point}" -z "${t_sample_ratio}" -s "${step_size}" -a "${smoothed_control_points}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}"

# "${convert_root_to_vtk}" -f "${smoothed_degenerate_point}" -o "${smoothed_degenerate_point}.csv" -i "${smoothed_mfa_file}" -d 0 -t 0

# "${tracking}" -f "${smoothed_mfa_file}" -b "${smoothed_tracking_result}" -z "${t_sample_ratio}" -g "${step_size}"  -a "${smoothed_control_points}" -x "${root_finding_epsilon}" -s "${smoothed_degenerate_point}" -p "${point_itr_threshold}"


# "${write_vtk}" -f "${smoothed_mfa_file}" -t "${smoothed_mfa_file}.vtk" -g 0 -z 0 -i "${data_type}" -u "${upsample_ratio}" #-s "${block}"

# python ./src/critical_point_tracking/time_data_convert.py -i "${smoothed_mfa_file}.vtk" -o "${smoothed_mfa_file}.vti"
# pvpython ./src/critical_point_tracking/extract_all_critical_points.py -i "${smoothed_mfa_file}.vti" -o "${smoothed_ttk_critical_point_file}.csv"

# pvpython ./src/critical_point_tracking/tracking_script.py -i "${mfa_file}.vti" -o "${ttk_tracking_file}"



##########################################################
# directly work on raw explicit function

for step_size in 4
do
#     "${degenerate_case_explicit}" -f "${data_type}" -b "${degenerate_point_original}_${step_size}.dat" -z "${step_size}" -s "${step_size}" -j "${J_threshold}" -p "${point_itr_threshold}" -g "${root_finding_epsilon}"

#     "${convert_root_to_vtk}" -f "${degenerate_point_original}_${step_size}.dat" -o "${degenerate_point_original}.csv" -j 0 -t 0

# #     # gdb --args 
    # "${tracking_explicit}" -f "${data_type}" -b "${ori_tracking_result}_${step_size}" -z "${step_size}" -g "${step_size}"  -x "${root_finding_epsilon}" -s "${degenerate_point_original}_${step_size}.dat" -p "${point_itr_threshold}" -e "${edge_type_file}_${step_size}.csv"


source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

    # python ./src/python/create_cube.py


source ~/enter/etc/profile.d/conda.sh
conda activate mfa

    # pvpython ./src/python/merge_obj_edge_type.py -i "${ori_tracking_result}_${step_size}.obj" -j "${edge_type_file}_${step_size}.csv" -o "${ori_tracking_result_w_type}_${step_size}.vtp"

# pvpython ./src/python/point_at_certain_time.py --input "${ori_tracking_result_w_type}_${step_size}.vtp" --output "${ori_tracking_result_w_type}_${step_size}_2.76.vtp" --z 2.76 --type-array EdgeValues

#     pvpython "${count_betti_num}" "${ori_tracking_result_w_type}_${step_size}.vtp"
done


#  python ./src/critical_point_tracking/time_data_convert.py -i "${ori_raw_data}.vtk" -o "${ori_raw_data}.vti"