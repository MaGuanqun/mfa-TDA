import pandas as pd
import argparse

def filter_points(input_file):
    # Read the CSV file
    points = pd.read_csv(input_file)
    
    # Filter points where z is 0
    filtered_points = points[points['x2'] == 0]
    
    return filtered_points

if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser(description='TTK-critical points.')
    parser.add_argument('-i', '--file_name', type=str, default='file_name.csv')

    args = parser.parse_args()
    
    input_file = args.file_name
    result = filter_points(input_file)
    
    print("filter points size ",len(result))
    
    result.to_csv(input_file, index=False)
    
    
