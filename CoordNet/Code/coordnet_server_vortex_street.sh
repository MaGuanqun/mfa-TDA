#!/bin/bash
#SBATCH -p spartacus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH -J coordnet
#SBATCH -o logs/%x-%j.out     # stdout file
#SBATCH -e logs/%x-%j.err     # stderr file


mkdir -p logs

export PYTHONUNBUFFERED=1

# source ~/.bashrc


application='super-spatial-temporal'
num_res=5
activate='sine' # sine, tanh
init_feature=64

num_epoch=300
omega=30.0


function_name=vortex_street_3d #vortex_street_3d, boussinesq_3d, fluid #cylinder2

raw_data="../Data/$function_name.bin"
raw_vti_file="../Result/$function_name/raw_file.vti"
raw_critical_point="../Result/$function_name/raw_critical_points.csv"


out_root=logs/$function_name
checkpoint=$out_root/checkpoints/model_final.pth


func_raw_data="../Result/$function_name/$application-$init_feature-$num_res"
vti_file="../Result/$function_name/$application-$init_feature-$num_res"
vti_critical_point="../Result/$function_name/vti_critical_points-$num_res"

convert_obj_domain=../../src/critical_point_tracking/convert_obj_back_to_domain.py



# python data_preprocessing.py
# 
python main.py --train 'train' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --lap_weight 0.0 --batch_size 16000 #--resume_epoch 270

# source ~/.bashrc

# source ~/enter/etc/profile.d/conda.sh
# conda activate siren
for step_size in 2 4 8 16 32 64
do
    python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000 --up_sample_ratio $step_size
done
# source ~/.bashrc
source /home/u1435513-gma/enter/etc/profile.d/conda.sh
conda activate mfa
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$CONDA_PREFIX/lib"

for step_size in 2 4 8 16 32 64
do
    python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}-${step_size}.dat" -o "${vti_file}-${step_size}.vti" -f "${function_name}" --step_size $step_size
    pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}-${step_size}.vti" -o "${vti_critical_point}-${step_size}.csv" -s 1
    # python $convert_obj_domain --csv "${vti_critical_point}-${step_size}.csv" --function "${function_name}" --back_to_ori 1
done




# python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000

# source ~/enter/etc/profile.d/conda.sh
# conda activate mfa_env


# source /home/u1435513-gma/enter/etc/profile.d/conda.sh
# conda activate mfa
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$CONDA_PREFIX/lib"
# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${func_raw_data}" -o "${vti_file}" -f "${function_name}" 
# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${vti_file}" -o "${vti_critical_point}"

# python ../../src/critical_point_tracking/binary_time_data_convert.py -i "${raw_data}" -o "${raw_vti_file}" -f "${function_name}"
# pvpython ../../src/critical_point_tracking/extract_all_critical_points.py -i "${raw_vti_file}" -o "${raw_critical_point}" -s 1


# tensorboard --logdir=logs/expotential/summaries/ --host=0.0.0.0
# ip addr show eth0 | grep "inet\b" | awk '{print $2}' | cut -d/ -f1
## after get ip address, open browser and go to http://ip_address:6006


# function_name=boussinesq_3d #quartic_potential_2, vortex_street,vortex_street_3d, hurricane_isabel,boussinesq_3d

# raw_data="../Data/$function_name.bin"
# raw_vti_file="../Result/$function_name/raw_file.vti"
# raw_critical_point="../Result/$function_name/raw_critical_points.csv"


# out_root=logs/$function_name
# checkpoint=$out_root/checkpoints/model_final.pth


# func_raw_data="../Result/$function_name/$application-$init_feature-$num_res.dat"
# vti_file="../Result/$function_name/$application-$init_feature-$num_res.vti"
# vti_critical_point="../Result/$function_name/vti_critical_points.csv"


# # python data_preprocessing.py
# # 
# python main.py --train 'train' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --lap_weight 0.0 --batch_size 16000 #--resume_epoch 270

# python main.py --train 'inf' --dataset $function_name --application $application --factor 1 --omega_0 $omega --init $init_feature --num_res $num_res --active $activate --num_epochs $num_epoch --batch_size 60000
