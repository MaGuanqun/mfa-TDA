source ~/enter/etc/profile.d/conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate siren




function_name=vortex_street_3d #quartic_potential_2, vortex_street,vortex_street_3d


omega=30.0

raw_data="../Data/$function_name.bin"
raw_vti_file="../Result/$function_name/raw_file.vti"
raw_critical_point="../Result/$function_name/raw_critical_points.csv"

application='super-spatial-temporal'
num_res=1
activate='sine' # sine, tanh
init_feature=64

num_epoch=300
out_root=logs/$function_name
checkpoint=$out_root/checkpoints/model_final.pth


func_raw_data="../Result/$function_name/$application-$init_feature-$num_res.dat"
vti_file="../Result/$function_name/$application-$init_feature-$num_res.vti"
vti_critical_point="../Result/$function_name/vti_critical_points.csv"


python main.py --train 'train' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --lap_weight 0.1

python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000

source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}" -o "${vti_file}" -f "${function_name}" 
# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}" -o "${vti_critical_point}"

# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${raw_data}" -o "${raw_vti_file}" -f "${function_name}"
# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${raw_vti_file}" -o "${raw_critical_point}"


# tensorboard --logdir=logs/expotential/summaries/ --host=0.0.0.0
# ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
## after get ip address, open browser and go to http://ip_address:6006