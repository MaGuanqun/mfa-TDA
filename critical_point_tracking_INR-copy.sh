#!/bin/bash
echo "Start Runing Script"


write_vtk="./build/src/convert/write_vtk"

tracking_INR="./build/src/critical_point_tracking/critical_point_tracking_INR"

degenerate_case_INR="./build/src/critical_point_tracking/degenerate_case_INR"
test_derivatives="./build/src/critical_point_tracking/test_derivatives"

convert_root_to_vtk="./build/src/critical_point/convert_root_to_vtk"


export_raw_data="./build/src/encode/analytical/export_raw_data"


control_point_smoothing="./build/src/critical_point_tracking/control_point_smoothing"

count_betti_num="./src/contour/count_betti_num.py"
convert_obj_domain=./src/critical_point_tracking/convert_obj_back_to_domain.py


step_size="32"
t_sample_ratio="32"
root_finding_epsilon="1e-2"
J_threshold="1e-1"

for data_type in "vortex_street_3d"
do
    echo "Processing data type: ${data_type}"
    # data_type="vortex_street_3d" #boussinesq_3d
    # data_type="quartic_potential_2"
    # data_type="sinc"
    save_folder="${data_type}_INR"

    control_points="./build/src/${save_folder}/${data_type}_cpt.dat"

    smoothed_control_points="./build/src/${save_folder}/${data_type}_cpt_smoothed.dat"


    mfa_file="./build/src/${save_folder}/${data_type}.mfa"

    smoothed_mfa_file="./build/src/${save_folder}/${data_type}_smoothed.mfa"

    ori_raw_data="./build/src/${save_folder}/${data_type}_raw.dat"

    degenerate_point_INR="./build/src/${save_folder}/${data_type}_degenerate_INR"

    boundary_point_INR="./build/src/${save_folder}/${data_type}_boundary"

    smoothed_degenerate_point="./build/src/${save_folder}/${data_type}_degenerate_smoothed.dat"

    tracking_result="./build/src/${save_folder}/${data_type}"
    smoothed_tracking_result="./build/src/${save_folder}/${data_type}_smoothed.obj"
    tracking_result_w_type="./build/src/${save_folder}/${data_type}"

    gradient_file="./build/src/${save_folder}/${data_type}.gradient.csv"

    edge_type_file="./build/src/${save_folder}/${data_type}_edge_type"

    ttk_tracking_file="./build/src/${save_folder}/ttk_${data_type}.vtu"
    ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt"

    smoothed_ttk_critical_point_file="./build/src/${save_folder}/ttk_${data_type}_cpt_smoothed.csv"


    # upsample_ratio="${step_size}-${step_size}-${t_sample_ratio}"

    #the reshold should be really small to raw explicit function

    point_itr_threshold="4.0"

    if [ "${data_type}" = "quartic_potential_2" ]; then
        input_model="./GridINR/SavedModels/${data_type}/fVSRN.pt"
    elif [ "${data_type}" = "vortex_street" ]; then
        input_model="./CoordNet/Exp/${data_type}/super-spatial-temporal-64-5.pt"
    elif [ "${data_type}" = "vortex_street_3d" ]; then
        input_model="./CoordNet/Exp/${data_type}/super-spatial-temporal-64-5.pt"
    elif [ "${data_type}" = "boussinesq_3d" ]; then
        input_model="./CoordNet/Exp/${data_type}/super-spatial-temporal-64-5.pt"
    fi


    step_size="32"

    # if [ "${data_type}" = "vortex_street_3d" ]; then
    #     "${tracking_INR}" -f "${data_type}" -b "${tracking_result}_${step_size}.obj" -z "${step_size}" -g "${step_size}"  -x "${root_finding_epsilon}" -s "${degenerate_point_INR}${step_size}.dat" -p "${point_itr_threshold}" -i "${input_model}" -e "${edge_type_file}_${step_size}.csv" -c 0 -j "${boundary_point_INR}"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     "${degenerate_case_INR}" -f "${data_type}" -b "${degenerate_point_INR}" -z "${step_size}" -s "${step_size}" -j "${J_threshold}" -g "${root_finding_epsilon}" -m "${input_model}"
    #     "${tracking_INR}" -f "${data_type}" -b "${tracking_result}_${step_size}.obj" -z "${step_size}" -g "${step_size}"  -x "${root_finding_epsilon}" -s "${degenerate_point_INR}${step_size}.dat" -p "${point_itr_threshold}" -i "${input_model}" -e "${edge_type_file}_${step_size}.csv" -c 1 -j "${boundary_point_INR}"
    # fi

    # 

    

    for step_size in "2" "4" "8" "16" "32"
    do
        # "${convert_root_to_vtk}" -f "${degenerate_point_INR}${step_size}.dat" -o "${degenerate_point_INR}_${step_size}.csv" -j 0 -t 1 
    #     echo "smoothing with step size: ${step_size}"
    #     t_sample_ratio="${step_size}"
        # "${tracking_INR}" -f "${data_type}" -b "${tracking_result}_${step_size}.obj" -z "${step_size}" -g "${step_size}"  -x "${root_finding_epsilon}" -s "${degenerate_point_INR}${step_size}.dat" -p "${point_itr_threshold}" -i "${input_model}" -e "${edge_type_file}_${step_size}.csv" -c 0 -j "${boundary_point_INR}"
        source ~/enter/etc/profile.d/conda.sh
        conda activate mfa
    #     python $convert_obj_domain --obj "${tracking_result}_${step_size}.obj" --function "${data_type}"
    #     python ./src/python/merge_obj_edge_type.py -i "${tracking_result}_${step_size}_rescale.obj" -j "${edge_type_file}_${step_size}.csv" -o "${tracking_result_w_type}_${step_size}.vtp"
        python "${count_betti_num}" "${tracking_result}_${step_size}.obj"
    done
done


if [ "${step_size}" = "16" ]; then
    # if [ "${data_type}" = "vortex_street_3d" ]; then
        spatial_step_size="0.00625" 
        t_step_size="0.00625"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     spatial_step_size="0.00694444"
    #     t_step_size="0.0110294"
    # fi

elif [ "${step_size}" = "8" ]; then
    # if [ "${data_type}" = "vortex_street_3d" ]; then
        spatial_step_size="0.0125" 
        t_step_size="0.0125"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     spatial_step_size="0.01388888"
    #     t_step_size="0.0220588"
    # fi
elif [ "${step_size}" = "4" ]; then
    # if [ "${data_type}" = "vortex_street_3d" ]; then
        spatial_step_size="0.025" 
        t_step_size="0.025"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     spatial_step_size="0.02777777"
    #     t_step_size="0.0441176"
    # fi
elif [ "${step_size}" = "32" ]; then
    # if [ "${data_type}" = "vortex_street_3d" ]; then
        spatial_step_size="0.003125" 
        t_step_size="0.003125"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     spatial_step_size="0.00347222"
    #     t_step_size="0.0055147"
    # fi
elif [ "${step_size}" = "2" ]; then
    # if [ "${data_type}" = "vortex_street_3d" ]; then
        spatial_step_size="0.05" 
        t_step_size="0.05"
    # elif [ "${data_type}" = "boussinesq_3d" ]; then
    #     spatial_step_size="0.00173611"
    #     t_step_size="0.00275735"
    # fi
fi



#2*step_size
# if [ "${step_size}" = "32" ]; then
#     # if [ "${data_type}" = "vortex_street_3d" ]; then
#         spatial_step_size="0.00625" 
#         t_step_size="0.00625"
#     # elif [ "${data_type}" = "boussinesq_3d" ]; then
#     #     spatial_step_size="0.00694444"
#     #     t_step_size="0.0110294"
#     # fi

# elif [ "${step_size}" = "16" ]; then
#     # if [ "${data_type}" = "vortex_street_3d" ]; then
#         spatial_step_size="0.0125" 
#         t_step_size="0.0125"
#     # elif [ "${data_type}" = "boussinesq_3d" ]; then
#     #     spatial_step_size="0.01388888"
#     #     t_step_size="0.0220588"
#     # fi
# elif [ "${step_size}" = "8" ]; then
#     # if [ "${data_type}" = "vortex_street_3d" ]; then
#         spatial_step_size="0.025" 
#         t_step_size="0.025"
#     # elif [ "${data_type}" = "boussinesq_3d" ]; then
#     #     spatial_step_size="0.02777777"
#     #     t_step_size="0.0441176"
#     # fi
# elif [ "${step_size}" = "4" ]; then
#     # if [ "${data_type}" = "vortex_street_3d" ]; then
#         spatial_step_size="0.05" 
#         t_step_size="0.05"
#     # elif [ "${data_type}" = "boussinesq_3d" ]; then
#     #     spatial_step_size="0.00173611"
#     #     t_step_size="0.00275735"
#     # fi
# elif [ "${step_size}" = "2" ]; then
#         spatial_step_size="0.1" 
#         t_step_size="0.1"
# fi


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


#gdb --args 

