source ~/enter/etc/profile.d/conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate siren




function_name=quartic_potential_2 #quartic_potential_2, vortex_street
output="./dataset/${function_name}.raw"


python generate_data.py --function $function_name --output $output 
