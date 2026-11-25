"""Join the two native h=0.125 crops to the h=0.5 Isabel candidate-2 boundary."""
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


FACE = np.dtype([('n', 'u1'), ('vertices', '<i4', (3,)), ('patch', 'u1')])


def read_ply(path):
    with path.open('rb') as stream:
        lines = []
        for line in stream:
            lines.append(line)
            if line == b'end_header\n':
                break
        if lines[:2] != [b'ply\n', b'format binary_little_endian 1.0\n'] or lines[-1] != b'end_header\n':
            raise ValueError(f'Expected a native binary boundary PLY in {path}')
        def count(element):
            return int(next(line for line in lines if line.startswith(b'element ' + element + b' ')).split()[-1])
        n, m = count(b'vertex'), count(b'face')
        meta = next(line for line in lines if line.startswith(b'comment fiber_volume_boundary ')).split()
        xyz = np.fromfile(stream, '<f8', 3*n).reshape(n, 3)
        cells = np.fromfile(stream, FACE, m)
        if len(cells) != m or not np.all(cells['n'] == 3):
            raise ValueError(f'Invalid faces in {path}')
    return xyz, cells['vertices'].copy(), cells['patch'].copy(), tuple(map(int, meta[3:6]))


def write_ply(path, xyz, faces, labels, corners):
    used, local = np.unique(faces, return_inverse=True)
    cells = np.empty(len(faces), FACE)
    cells['n'], cells['vertices'], cells['patch'] = 3, local.reshape(-1, 3), labels
    with path.open('wb') as stream:
        metadata = ' '.join(map(str, corners))
        stream.write((f'ply\nformat binary_little_endian 1.0\ncomment fiber_volume_boundary 0 {metadata} 0\n'
                      f'element vertex {len(used)}\nproperty double x\nproperty double y\nproperty double z\n'
                      f'element face {len(faces)}\nproperty list uchar int vertex_indices\n'
                      'property uchar patch\nend_header\n').encode())
        xyz[used].astype('<f8').tofile(stream)
        cells.tofile(stream)


def cut_loops(faces, labels, retained, donor=False, expected=(1, 2)):
    """Find the oriented border of the retained faces, with one loop per side."""
    edges = np.vstack([faces[:, [i, (i+1) % 3]] for i in range(3)])
    keys = edges.min(axis=1).astype('i8') * (int(faces.max())+1) + edges.max(axis=1)
    order = np.argsort(keys)
    _, first, count = np.unique(keys[order], return_index=True, return_counts=True)
    if np.any(count > 2):
        raise ValueError('Cut contains a nonmanifold edge')
    pair = first[count == 2]
    a, b = order[pair], order[pair+1]
    use_a, use_b = retained[a % len(faces)], retained[b % len(faces)]
    border = np.where(use_a, a, b)[use_a != use_b]
    if donor:
        singles = order[first[count == 1]]
        if np.any(retained[singles % len(faces)]):
            raise ValueError('Refined patch extends to an unpaired donor edge')
    successors = {u: (v, label) for (u, v), label in zip(edges[border], labels[border % len(faces)])}
    if len(successors) != len(border) or set(successors) != set(edges[border, 1]):
        raise ValueError('Cut produced a branched or unoriented rim')
    loops = {}
    while successors:
        start = v = next(iter(successors))
        vertices, sides = [], set()
        while v in successors:
            vertices.append(v)
            v, side = successors.pop(v)
            sides.add(side)
        if v != start or len(sides) != 1 or side in loops:
            raise ValueError('Expected a single closed loop on each range side')
        loops[side] = np.array(vertices)
    if expected is not None and set(loops) != set(expected):
        raise ValueError('Expected one rim on each requested side')
    return loops


def join_strip(xyz, outer, inner, ref_points, ref_normals, step):
    # A step in this lattice advances by one edge of one rim and adds a triangle.
    # Minimize the total cost over consistently oriented paths around both rims.
    outer = outer[::-1]
    d = np.linalg.norm(xyz[outer, None] - xyz[inner], axis=2)
    p, q = np.unravel_index(np.argmin(d), d.shape)
    outer, inner = np.roll(outer, -p), np.roll(inner, -q)
    na, nb = len(outer), len(inner)
    a = np.r_[outer, outer[0]]
    tree = cKDTree(ref_points)

    def cost(triangles):
        coords = xyz[triangles]
        u, v = coords[..., 1, :] - coords[..., 0, :], coords[..., 2, :] - coords[..., 0, :]
        cross = np.cross(u, v)
        area2 = np.linalg.norm(cross, axis=-1)
        centroid = coords.mean(axis=-2)
        nearest = tree.query(centroid.reshape(-1, 3))[1]
        cosine = (cross * ref_normals[nearest].reshape(cross.shape)).sum(axis=-1) / np.maximum(area2, 1e-20)
        longest = np.maximum.reduce([np.linalg.norm(e, axis=-1) for e in (u, v, u-v)])
        score = longest**2 / np.maximum(area2, 1e-20) + 4*(1-cosine) + longest**2
        return np.where((cosine > .10) & (longest <= 2*step) & (area2 > 1e-10), score, np.inf)

    best = None
    for shift in [0, -1, 1, -2, 2]:
        b = np.roll(inner, shift)
        b = np.r_[b, b[0]]
        i, j = np.meshgrid(np.arange(na), np.arange(nb+1), indexing='ij')
        advance_outer = np.stack([a[i], a[i+1], b[j]], axis=-1)
        i, j = np.meshgrid(np.arange(na+1), np.arange(nb), indexing='ij')
        advance_inner = np.stack([a[i], b[j+1], b[j]], axis=-1)
        co, ci = cost(advance_outer), cost(advance_inner)
        scores = np.full((na+1, nb+1), np.inf)
        scores[0, 0] = 0
        from_outer = np.zeros(scores.shape, bool)
        for i in range(na+1):
            for j in range(nb+1):
                if i+j == 0:
                    continue
                u = scores[i-1, j] + co[i-1, j] if i else np.inf
                v = scores[i, j-1] + ci[i, j-1] if j else np.inf
                from_outer[i, j], scores[i, j] = u < v, min(u, v)
        score = scores[-1, -1]
        if not np.isfinite(score) or (best is not None and best[0] <= score):
            continue
        i, j = na, nb
        faces = []
        while i or j:
            if from_outer[i, j]:
                i -= 1
                faces.append(advance_outer[i, j])
            else:
                j -= 1
                faces.append(advance_inner[i, j])
        best = score, np.array(faces)
    if best is None:
        raise ValueError('No oriented joining strip within the edge-length limit')
    return best[1]


def join(coarse, isolated, main, prefix):
    x, f, labels, corners = read_ply(coarse)
    y, shell, shell_labels, _ = read_ply(isolated)
    low, high = np.array([323.5, 351.5, 42.25]), np.array([326.5, 354.5, 43.25])
    coords = x[f]
    inside = np.all((coords >= low) & (coords <= high), axis=(1, 2))
    overlap = np.all((coords.max(axis=1) >= low) & (coords.min(axis=1) <= high), axis=1)
    if (not inside.any() or not np.array_equal(inside, overlap) or
            not np.all((y > low) & (y < high)) or np.any(shell_labels > 3)):
        raise ValueError('The small shell must be isolated from the crop boundary')
    f = np.vstack([f[~inside], shell + len(x)])
    x = np.vstack([x, y])
    labels = np.r_[labels[~inside], shell_labels]

    y, fine, fine_labels, _ = read_ply(main)
    offset = len(x)
    xyz = np.vstack([x, y])
    coarse_centers, fine_centers = x[f].mean(axis=1), y[fine].mean(axis=1)
    center = np.array([279.1, 350.1])
    r = np.linalg.norm(coarse_centers[:, :2] - center, axis=1)
    rf = np.linalg.norm(fine_centers[:, :2] - center, axis=1)
    nearby = (r < 3) & (coarse_centers[:, 2] > 41.8) & (coarse_centers[:, 2] < 44)
    remove = nearby & (((labels == 1) & (r < 1.65)) |
                       ((labels == 2) & (r < 2.5) & (coarse_centers[:, 2] > 42.5)))
    take = ((fine_labels == 1) & (rf < 1.25)) | ((fine_labels == 2) & (fine_centers[:, 2] > 42.8))
    outer = cut_loops(f[nearby], labels[nearby], ~remove[nearby])
    inner = cut_loops(fine+offset, fine_labels, take, donor=True)
    strips, side_ids = [], []
    for label in (1, 2):
        ref = fine_labels == label
        p = y[fine[ref]]
        normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        strip = join_strip(xyz, outer[label], inner[label], fine_centers[ref], normal, .5)
        strips.append(strip)
        side_ids.append(np.full(len(strip), label))
    strips, side_ids = np.vstack(strips), np.concatenate(side_ids)
    write_ply(prefix.with_suffix('.ply'), xyz, np.vstack([f[~remove], fine[take]+offset, strips]),
              np.r_[labels[~remove], fine_labels[take], side_ids], corners)
    write_ply(prefix.with_name(prefix.name+'_transition.ply'), xyz, strips, side_ids, (0, 0, 0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['coarse', 'isolated', 'main', 'prefix']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    join(args.coarse, args.isolated, args.main, args.prefix)
