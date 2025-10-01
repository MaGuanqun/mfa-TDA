source ~/enter/etc/profile.d/conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate siren


activation_type=tanh
num_epoch=1000
hidden_feature=64
hidden_layer=16
function_name=schwefel #expotential
model_mode=mlp # mlp,nerf

python experiment_scripts/train_function.py --model_type=$activation_type --experiment_name=expotential --function=$function_name --num_epochs=$num_epoch --hidden_features=$hidden_feature --num_hidden_layers=$hidden_layer --mode=$model_mode
# python experiment_scripts/train_img.py --model_type=$activation_type --experiment_name=expotential

python experiment_scripts/compute_h.py --hidden_features=$hidden_feature --num_hidden_layers=$hidden_layer --model_type=$activation_type --function=$function_name --mode=$model_mode


# tensorboard --logdir=logs/expotential/summaries/ --host=0.0.0.0
# ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
## after get ip address, open browser and go to http://ip_address:6006