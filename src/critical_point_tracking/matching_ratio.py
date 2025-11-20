import argparse
import numpy as np

def load_obj_vertices(obj_path):
    """Load vertex positions (x, y, z) from an .obj file."""
    vertices = []
    with open(obj_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                parts = line.split()
                if len(parts) >= 4:
                    x, y, z = map(float, parts[1:4])
                    vertices.append((x, y, z))
    return np.array(vertices, dtype=float)


def load_csv_points(csv_path, has_header=True):
    """Load (x, y, z) points from a CSV file."""
    skip = 1 if has_header else 0
    return np.loadtxt(csv_path, delimiter=',', skiprows=skip, usecols=(0, 1, 2))


def compute_match_ratio(obj_pts, csv_pts, xy_threshold, z_threshold, csv_has_header):
    if obj_pts.size == 0:
        raise ValueError("No vertices found in the OBJ file.")
    if csv_pts.size == 0:
        raise ValueError("No points found in the CSV file.")

    # Try KD-tree acceleration
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(obj_pts[:, :2])
        use_kdtree = True
        print("Using KD-tree for nearest neighbor search.")
    except ImportError:
        use_kdtree = False
        print("scipy not available; using slower method for nearest neighbor search.")

    matched = 0
    n_csv = csv_pts.shape[0]

    if use_kdtree:
        neighbors_list = tree.query_ball_point(csv_pts[:, :2], r=xy_threshold)

        for i, neighbors in enumerate(neighbors_list):
            if not neighbors:
                continue
            z_csv = csv_pts[i, 2]
            if np.any(np.abs(obj_pts[neighbors, 2] - z_csv) < z_threshold):
                matched += 1
    else:
        # Slower fallback
        obj_xy = obj_pts[:, :2]
        obj_z = obj_pts[:, 2]

        for x, y, z in csv_pts:
            diff_xy = obj_xy - np.array([x, y])
            dist_xy = np.sqrt(np.sum(diff_xy**2, axis=1))
            candidates = dist_xy < xy_threshold
            if np.any(np.abs(obj_z[candidates] - z) < z_threshold):
                matched += 1

    ratio = matched / n_csv
    print(f"Matched: {matched} / {n_csv}")
    print(f"Ratio:   {ratio:.6f}")

    return ratio


def compute_match_ratio_different(obj_path, csv_path, xy_threshold, z_threshold, csv_has_header):
    obj_pts = load_obj_vertices(obj_path)
    csv_pts = load_csv_points(csv_path, has_header=csv_has_header)
    print("discrete result on our result:")
    compute_match_ratio(obj_pts, csv_pts, xy_threshold, z_threshold, csv_has_header)
    print("our result on discrete result:")
    compute_match_ratio(csv_pts, obj_pts, xy_threshold, z_threshold, csv_has_header)

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare OBJ vertices and CSV points.")
    parser.add_argument("obj_file", type=str, help="Path to OBJ file.")
    parser.add_argument("csv_file", type=str, help="Path to CSV file.")
    parser.add_argument("xy_threshold", type=float, help="XY-plane distance threshold.")
    parser.add_argument("z_threshold", type=float, help="Z distance threshold.")
    parser.add_argument("--header", action="store_true",
                        help="Indicate that CSV has a header row.")

    args = parser.parse_args()

    compute_match_ratio_different(
        obj_path=args.obj_file,
        csv_path=args.csv_file,
        xy_threshold=args.xy_threshold,
        z_threshold=args.z_threshold,
        csv_has_header=args.header,
    )