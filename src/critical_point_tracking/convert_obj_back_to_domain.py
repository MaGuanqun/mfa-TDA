import numpy as np
import csv
import argparse
import os

def new_domain(function_name='vortex_street_3d'):

    if function_name == 'vortex_street_3d':
        return np.array([0.0, 639.0, 0.0, 79.0, 1350.0, 1500.0])
    elif function_name == 'boussinesq_3d':
        return np.array([0.0,149.0,0.0,449.0,0.0,199.0])
    elif function_name == 'fluid':
        return np.array([0.0,511.0,0.0,511.0,0.0,99.0])
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")
    
def ori_domain(function_name='vortex_street_3d'):

    if function_name == 'vortex_street_3d':
        return np.array([-0.5, 7.5, -0.5, 0.5, 13.5,15.0])
    elif function_name == 'boussinesq_3d':
        return np.array([-0.5,0.5,-0.5,2.5,0.0,2.0])
    elif function_name == 'fluid':
        return np.array([0,1,0,1,0,1])
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")
    
def convert_point_to_new_domain(input_points, function_name='vortex_street_3d'):
    """
    Map points from ori_domain -> new_domain.
    input_points: (N, 3) array in original domain coordinates.
    """
    ori_dom = ori_domain(function_name)
    new_dom = new_domain(function_name)

    # Scale = (new_range / ori_range)
    scale_x = (new_dom[1] - new_dom[0]) / (ori_dom[1] - ori_dom[0])
    scale_y = (new_dom[3] - new_dom[2]) / (ori_dom[3] - ori_dom[2])
    scale_z = (new_dom[5] - new_dom[4]) / (ori_dom[5] - ori_dom[4])

    converted_points = np.zeros_like(input_points, dtype=float)
    converted_points[:, 0] = new_dom[0] + (input_points[:, 0] - ori_dom[0]) * scale_x
    converted_points[:, 1] = new_dom[2] + (input_points[:, 1] - ori_dom[2]) * scale_y
    converted_points[:, 2] = new_dom[4] + (input_points[:, 2] - ori_dom[4]) * scale_z

    return converted_points

def convert_back_to_ori_domain(input_points, function_name='vortex_street_3d'):
    """
    Map points from new_domain -> ori_domain.
    input_points: (N, 3) array in new domain coordinates.
    """
    ori_dom = ori_domain(function_name)
    new_dom = new_domain(function_name)

    # Scale = (ori_range / new_range)
    scale_x = (ori_dom[1] - ori_dom[0]) / (new_dom[1] - new_dom[0])
    scale_y = (ori_dom[3] - ori_dom[2]) / (new_dom[3] - new_dom[2])
    scale_z = (ori_dom[5] - ori_dom[4]) / (new_dom[5] - new_dom[4])

    converted_points = np.zeros_like(input_points, dtype=float)
    converted_points[:, 0] = ori_dom[0] + (input_points[:, 0] - new_dom[0]) * scale_x
    converted_points[:, 1] = ori_dom[2] + (input_points[:, 1] - new_dom[2]) * scale_y
    converted_points[:, 2] = ori_dom[4] + (input_points[:, 2] - new_dom[4]) * scale_z

    return converted_points


def remap_obj_vertices(input_obj_path, output_obj_path, function_name='vortex_street_3d', back_to_ori=0):
    """Read OBJ, map all vertex positions back to ori_domain, and write new OBJ."""
    with open(input_obj_path, 'r') as f:
        lines = f.readlines()

    vertex_positions = []
    vertex_line_indices = []

    # 1) Collect all 'v' vertex lines
    for i, line in enumerate(lines):
        # Standard vertex lines start with "v " (not "vn", "vt", etc.)
        if line.startswith('v '):
            parts = line.strip().split()
            if len(parts) >= 4:
                x, y, z = map(float, parts[1:4])
                vertex_positions.append([x, y, z])
                vertex_line_indices.append(i)

    if not vertex_positions:
        print("No vertex positions ('v' lines) found in OBJ.")
        return

    vertex_positions = np.array(vertex_positions, dtype=float)

    # 2) Apply your domain conversion
    if back_to_ori==0:
        converted_positions = convert_point_to_new_domain(
            vertex_positions,
            function_name=function_name
        )
    else:
        converted_positions = convert_back_to_ori_domain(
            vertex_positions,
            function_name=function_name
        )

    # 3) Replace the vertex lines with converted coordinates
    for idx, (x, y, z) in zip(vertex_line_indices, converted_positions):
        parts = lines[idx].strip().split()
        # Keep any extra data on the vertex line (e.g. colors, w)
        extra = parts[4:]
        new_line = f"v {x:.6f} {y:.6f} {z:.6f}"
        if extra:
            new_line += " " + " ".join(extra)
        new_line += "\n"
        lines[idx] = new_line

    # 4) Write new OBJ
    with open(output_obj_path, 'w') as f:
        f.writelines(lines)

    print(f"Converted OBJ saved to: {output_obj_path}")



def remap_csv_points(input_csv_path, output_csv_path, function_name='vortex_street_3d', back_to_ori=0):
    """
    Read CSV of points (x,y,z), map ori_domain -> new_domain, write new CSV.
    Keeps the first header row.
    """
    rows = []
    header = None
    first_row = True

    with open(input_csv_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue

            if first_row:
                # save header exactly as-is
                header = row
                first_row = False
                continue

            # convert numeric rows only
            try:
                x, y, z = map(float, row[:3])
                rows.append([x, y, z])
            except ValueError:
                # If a non-numeric data row appears AFTER header → skip
                continue

    if not rows:
        print("[ERROR] CSV contains no numeric rows with x,y,z.")
        return

    data = np.array(rows, dtype=float)
    if back_to_ori==0:
        converted = convert_point_to_new_domain(data, function_name)
    else:
        converted = convert_back_to_ori_domain(data, function_name)

    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)        # ✓ write header back
        for x, y, z in converted:
            writer.writerow([f"{x:.6f}", f"{y:.6f}", f"{z:.6f}"])

    print(f"[OK] Converted CSV saved → {output_csv_path}")


    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Convert OBJ or CSV points back to original domain.")

    parser.add_argument("--obj", type=str, default=None, help="Input OBJ file")
    parser.add_argument("--csv", type=str, default=None, help="Input CSV file")
    parser.add_argument("--function", type=str, default="vortex_street_3d",
                        help="Function name for domain conversion")
    parser.add_argument("--back_to_ori", type=int, default=0,
                        help="If 1, convert from new_domain back to ori_domain")
    
    args = parser.parse_args()

    # ---------------------- OBJ Processing ----------------------
    if args.obj:
        base = os.path.splitext(args.obj)[0]
        default_output_obj = base + "_rescale.obj"
        remap_obj_vertices(args.obj, default_output_obj, args.function, args.back_to_ori)

    # ---------------------- CSV Processing ----------------------
    if args.csv:
        base = os.path.splitext(args.csv)[0]
        default_output_csv = base + "_rescale.csv"
        remap_csv_points(args.csv, default_output_csv, args.function, args.back_to_ori)