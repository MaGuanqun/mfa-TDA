source ~/enter/etc/profile.d/conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate siren


activation_type=sine # sine, tanh
num_epoch=3000
hidden_feature=512
hidden_layer=10
function_name=quartic_potential_2 #expotential #schwefel #quartic_potential_2
model_mode=mlp # mlp,nerf
omega=20.0

out_root=logs/$function_name
checkpoint=$out_root/checkpoints/model_final.pth


func_raw_data="logs/$function_name/f.bin"
vti_file="logs/$function_name/f.vti"
vti_critical_point="logs/$function_name/vti_critical_points.csv"

raw_data="logs/$function_name/f_true.bin"
raw_vti_file="logs/$function_name/f_true.vti"
raw_critical_point="logs/$function_name/raw_critical_points.csv"

# python experiment_scripts/train_function.py --model_type=$activation_type --experiment_name=$function_name --function=$function_name --num_epochs=$num_epoch --hidden_features=$hidden_feature --num_hidden_layers=$hidden_layer --mode=$model_mode --omega=$omega
# python experiment_scripts/train_img.py --model_type=$activation_type --experiment_name=expotential

python experiment_scripts/sample_model.py --hidden_features=$hidden_feature --num_hidden_layers=$hidden_layer --model_type=$activation_type --function=$function_name --mode=$model_mode --omega=$omega --checkpoint=$checkpoint --out_root=$out_root


source ~/enter/etc/profile.d/conda.sh
conda activate mfa_env

python ../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}" -o "${vti_file}" -f "${function_name}"
pvpython ../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}" -o "${vti_critical_point}"

# python ../src/critical_point_tracking/binary_time_data_convert.py -i "${raw_data}" -o "${raw_vti_file}" -f "${function_name}"
# pvpython ../src/critical_point_tracking/extract_all_critical_points.py -i "${raw_vti_file}" -o "${raw_critical_point}"


# tensorboard --logdir=logs/expotential/summaries/ --host=0.0.0.0
# ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
## after get ip address, open browser and go to http://ip_address:6006