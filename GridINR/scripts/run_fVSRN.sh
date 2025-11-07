#! /bin/bash


function_name=quartic_potential_2 #quartic_potential_2, vortex_street


if [ "$function_name" = "quartic_potential_2" ]; then
    raw_data="../generate_dataset/dataset/quartic_potential_2.raw"
    data_dim="100,100,100"
    feature_grid_shape="64,64,64"
elif [ "$function_name" = "vortex_street" ]; then
    raw_data="./Data/$function_name.bin"
    data_dim="100,80,50"
    feature_grid_shape="64,64,64"
fi


raw_vti_file="./Result/$function_name.raw_file.vti"
raw_critical_point="./Result/$function_name.raw_critical_points.csv"

func_raw_data="./Output/$function_name.raw"
vti_file="./Result/$function_name.vti"
vti_critical_point="./Result/${function_name}_critical_points.csv"



cd ..

# # training
# python3 train.py --model_type fVSRN --data_path $raw_data --dataset_name $function_name --data_dims $data_dim --num_positional_encoding_terms 0 --feature_grid_shape $feature_grid_shape --second_deriv_weight 2e-7
                 
# num_positional_encoding_terms: you may probably prefer 0 if you want a smooth reconstruction
# feature_grid_shape: increase this value to capture more details, but may not be smooth

# inference
python3 infer.py --model_type fVSRN --load_from ./SavedModels/$function_name


source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

python ../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}" -o "${vti_file}" -f "${function_name}" 
pvpython ../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}" -o "${vti_critical_point}"


# python ../src/critical_point_tracking/binary_time_data_convert.py -i "${raw_data}" -o "${raw_vti_file}" -f "${function_name}" --float_type float64
# pvpython ../src/critical_point_tracking/extract_all_critical_points.py -i "${raw_vti_file}" -o "${raw_critical_point}"