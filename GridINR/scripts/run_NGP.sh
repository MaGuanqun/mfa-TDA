#! /bin/bash
cd ..
# training
python3 train.py --model_type NGP \
                 --data_path ./Data/vortex_street.bin \
                 --dataset_name vortex_street \
                 --data_dims 100,80,50
# inference
python3 infer.py --load_from ./SavedModels/vortex_street


