#! /bin/bash
cd ..
# training
python3 train.py --model_type fVSRN \
                 --data_path ./Data/vortex_street.bin \
                 --dataset_name vortex_street \
                 --data_dims 100,80,50 \
                 --num_positional_encoding_terms 0 \
                 --feature_grid_shape 64,64,64 \
                 
# num_positional_encoding_terms: you may probably prefer 0 if you want a smooth reconstruction
# feature_grid_shape: increase this value to capture more details, but may not be smooth

# inference
python3 infer.py --model_type fVSRN \
                 --load_from ./SavedModels/vortex_street