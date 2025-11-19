source ~/enter/etc/profile.d/conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate siren




function_name=boussinesq_3d #vortex_street_3d, boussinesq_3d


omega=30.0

raw_data="../Data/${function_name}.bin"
raw_vti_file="../Result/$function_name/raw_file.vti"
new_raw_vtk_file="../Result/$function_name/new_raw_file.vtk"
new_raw_vti_file="../Result/$function_name/new_raw_file.vti"
new_raw_data="../Data/$function_name.bin"
raw_critical_point="../Result/$function_name/raw_critical_points.csv"
new_raw_critical_point="../Result/$function_name/new_raw_critical_points.csv"

application='super-spatial-temporal'
num_res=5
activate='sine' # sine, tanh
init_feature=64

num_epoch=300
out_root=logs/$function_name
checkpoint=$out_root/checkpoints/model_final.pth


func_raw_data="../Result/$function_name/$application-$init_feature-$num_res.dat"
ori_vti_file="../Result/$function_name/$application-$init_feature-$num_res.ori_vti"
vti_file="../Result/$function_name/$application-$init_feature-$num_res.vti"
vti_critical_point="../Result/$function_name/vti_critical_points.csv"


# python data_preprocessing.py

# python main.py --train 'train' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --lap_weight 0.0 --batch_size 16000

# python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000

# for step_size in 2 4 8 16 32
# do
    # python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000 --up_sample_ratio $step_size
# done

source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}" -o "${vti_file}" -f "${function_name}"



# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}" -o "${vti_critical_point}"




python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${raw_data}" -o "${raw_vti_file}" -f "${function_name}"
# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${raw_vti_file}" -o "${raw_critical_point}"
# # # conda deactivate

# pvpython ../../src/critical_point_tracking/persistence_simplification.py -i "${raw_vti_file}" -o "${new_raw_vti_file}" -b "${new_raw_data}" -s 2.0 --median_radius 0


# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${new_raw_data}" -o "${new_raw_vti_file}" -f "${function_name}"

# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${new_raw_vti_file}" -o "${new_raw_critical_point}"


# convert_obj_domain=../../src/critical_point_tracking/convert_obj_back_to_domain.py
# python $convert_obj_domain --csv "${new_raw_critical_point}" --function "${function_name}"


# tensorboard --logdir=logs/expotential/summaries/ --host=0.0.0.0
# ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
## after get ip address, open browser and go to http://ip_address:6006