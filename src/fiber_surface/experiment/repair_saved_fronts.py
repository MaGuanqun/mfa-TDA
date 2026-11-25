"""Repair a saved native mesh using its own vertices, faces and completed corner traces."""
import argparse
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from join_isabel_refinements import cut_loops, join_strip, read_ply, write_ply
from ttk_geometry import original_curves


def defects(faces, vertex_count):
    """Find open rims and vertices with disconnected incident triangle fans."""
    m = len(faces)
    if not m:
        return [], np.empty(0, dtype=int)
    corners = faces.T.ravel()
    successor = (np.arange(3*m) + m) % (3*m)
    other = corners[successor]
    edge = np.minimum(corners, other).astype('i8') * vertex_count + np.maximum(corners, other)
    order = np.argsort(edge)
    _, first, count = np.unique(edge[order], return_index=True, return_counts=True)
    if (count > 2).any():
        raise ValueError('Expected edges with at most two incident faces')
    singles = order[first[count == 1]]
    paired = first[count == 2]
    a, b = order[paired], order[paired+1]
    if not np.all(corners[a] == other[b]):
        raise ValueError('Inconsistent edge orientation')
    # Join the two face corners at each end of every interior edge. Their
    # connected components give exactly the triangle fans of each vertex.
    left, right = np.r_[a, successor[a]], np.r_[successor[b], b]
    graph = coo_matrix((np.ones(len(left), dtype='b'), (left, right)), shape=(3*m, 3*m))
    nfan, fan = connected_components(graph, directed=False)
    unique = np.unique(corners.astype('i8') * nfan + fan)
    pinches = np.flatnonzero(np.bincount(unique // nfan, minlength=vertex_count) > 1)
    boundary = np.column_stack([corners[singles], other[singles]])
    affected = np.union1d(boundary.ravel(), pinches)
    if not len(affected):
        return [], pinches
    edges = np.searchsorted(affected, boundary)
    graph = coo_matrix((np.ones(len(edges), dtype='b'), edges.T), shape=(len(affected), len(affected)))
    n, groups = connected_components(graph, directed=False)
    return [affected[groups == i] for i in range(n)], pinches


def unit(vectors):
    return vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)


def loop_orders(xyz, ids):
    yield ids, 'saved'
    points = xyz[ids] - xyz[ids].mean(axis=0)
    basis = np.linalg.svd(points, full_matrices=False)[2][:2]
    uv = points @ basis.T
    around = np.argsort(np.arctan2(uv[:, 1], uv[:, 0]))
    # Only introduce an alternative order if the saved polygon folds across
    # itself in its best-fitting plane. Keep each exact corner sample in either case.
    def crosses(a, b, c, d):
        system = np.stack((b-a, c-d)).T
        if abs(np.linalg.det(system)) <= 1e-12 * np.linalg.norm(b-a) * np.linalg.norm(d-c):
            return False
        s = np.linalg.solve(system, c-a)
        return np.all((s > 0) & (s < 1))
    n = len(ids)
    if any(crosses(uv[i], uv[(i+1) % n], uv[j], uv[(j+1) % n])
           for i in range(n) for j in range(i+2, n) if (j+1) % n != i):
        yield ids[around], 'untangled'


def audit(bin_dir, stages, path, prefix, *, exhaustive):
    command = [bin_dir/'fiber_grid_recovery', 'audit', stages, path, prefix]
    if exhaustive:
        command.append('--all-face-samples')
    with prefix.with_suffix('.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    return json.loads(prefix.with_suffix('.json').read_text())


def topology(report):
    return {key: int(value) for key, value in re.findall(r'(\w+)=(-?\d+)', report['topology'])}


def valid_local(result, boundary_edges, angle, tolerance):
    counts = topology(result)
    return (result['intersections_checked'] and counts['boundary_edges'] == boundary_edges
            and all(v == 0 for k, v in counts.items() if k not in ('boundary_edges', 'Euler'))
            and result['seam_failures'] == result['seam_chord_failures'] == result['face_sites_projection_failed'] == 0
            and result['face_sites_max_angle'] <= angle
            and result['face_sites_max_projection'] <= tolerance)


def disk_candidates(xyz, faces, labels, selected, ring, label, face_centers, face_normals):
    """Preserve the trace and try interior vertices already sampled on this side."""
    vertices = np.setdiff1d(np.unique(faces[selected & (labels == label)]), ring)
    tree = cKDTree(face_centers[labels == label])
    normals = face_normals[labels == label]
    choices = []
    for v in vertices:
        triangles = np.stack((ring, np.roll(ring, -1), np.full(len(ring), v)), axis=1)
        p = xyz[triangles]
        n = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])
        if (np.linalg.norm(n, axis=1) < 1e-14).any():
            continue
        cosine = (unit(n) * normals[tree.query(p.mean(axis=1))[1]]).sum(axis=1)
        if (cosine <= 0).any():
            continue
        choices.append((cosine.min(), triangles))
    return [t for _, t in sorted(choices, key=lambda row: -row[0])[:4]]


def proposals(xyz, faces, labels, selected, ring, sides, max_edge):
    outer = cut_loops(faces, labels, ~selected, expected=None)
    if outer and set(outer) != set(sides):
        return
    p = xyz[faces]
    centers = p.mean(axis=1)
    normals = unit(np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0]))
    for order, method in loop_orders(xyz, ring):
        for sign in (1, -1):
            alternatives = []
            for k, label in enumerate(sides):
                inner = order[::sign if k == 0 else -sign]
                if not outer:
                    alternatives.append(disk_candidates(xyz, faces, labels, selected, inner, label, centers, normals))
                    continue
                use = labels == label
                try:
                    strip = join_strip(xyz, outer[label], inner, centers[use], normals[use], max_edge/2)
                    alternatives.append([strip])
                except ValueError:
                    alternatives.append([])
            for pair in product(*alternatives):
                t = np.vstack(pair)
                longest = max(np.linalg.norm(xyz[t[:, j]]-xyz[t[:, (j+1) % 3]], axis=1).max()
                              for j in range(3))
                if longest <= max_edge:
                    yield t, np.repeat(sides, [len(side) for side in pair]), method, float(longest)


def fit_replacement(bin_dir, stages, work, xyz, local_faces, local_labels, t, sides,
                    corners, expected_edges, angle, tolerance, levels):
    """Keep exposed rim vertices fixed and accept whole conforming refinement levels."""
    if not levels:
        return xyz, t, sides, 0
    base = work/'replacement'
    write_ply(base.with_suffix('.ply'), xyz, t, sides, corners)
    previous = audit(bin_dir, stages, base.with_suffix('.ply'), base, exhaustive=True)
    originals = np.unique(t)
    source = base.with_suffix('.ply')
    chosen_xyz, chosen_t, chosen_sides, accepted = xyz, t, sides, 0
    for level in range(1, levels+1):
        target = work/f'fit_{level}.ply'
        with target.with_suffix('.log').open('w') as log:
            subprocess.run([bin_dir/'fiber_refine_saved_patch', stages, source, target],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        points, triangles, labels, _ = read_ply(target)
        if len(triangles) == len(chosen_t):
            break
        if not np.array_equal(points[:len(originals)], xyz[originals]):
            raise RuntimeError('Patch refinement changed an original vertex')
        ref = audit(bin_dir, stages, target, work/'fit_audit', exhaustive=True)
        if ref['projection_p95'] >= previous['projection_p95']:
            break
        candidate_x = np.vstack([xyz, points[len(originals):]])
        ids = np.r_[originals, np.arange(len(xyz),len(candidate_x))]
        candidate_t = ids[triangles]
        joined = work/'fit_with_neighbors'
        write_ply(joined.with_suffix('.ply'), candidate_x, np.vstack([local_faces,candidate_t]),
                  np.r_[local_labels,labels], corners)
        result = audit(bin_dir, stages, joined.with_suffix('.ply'), joined, exhaustive=True)
        if not valid_local(result, expected_edges, angle, tolerance):
            break
        chosen_xyz, chosen_t, chosen_sides, accepted = candidate_x, candidate_t, labels, level
        source, previous = target, ref
    return chosen_xyz, chosen_t, chosen_sides, accepted


def reconstruct(run, out, bin_dir, settings):
    stages = run/'stages'
    config = (stages/'run.conf').read_text().splitlines()
    opt = next(s.split()[1:] for s in config if s.startswith('options '))
    h, tolerance = float(opt[0]), float(opt[3])
    angle = settings.normal_angle or float(opt[7])
    bmin, bmax = [np.fromstring(next(s for s in config if s.startswith(key)).split(' ', 1)[1], sep=' ')
                  for key in ('box_min ', 'box_max ')]
    out.mkdir(parents=True, exist_ok=False)
    ordinary = out/'rim'
    with (out/'rim.log').open('w') as log:
        subprocess.run([bin_dir/'fiber_repair_boundary', stages, run/'volume_boundary.ply', ordinary, '--repair'],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    xyz, faces, labels, corners = read_ply(ordinary.with_suffix('.ply'))
    groups, pinches = defects(faces, len(xyz))
    before = audit(bin_dir, stages, ordinary.with_suffix('.ply'), out/'rim_audit', exhaustive=False)
    curves = [(k, i, path[:-1]) for k, loops in enumerate(original_curves(stages))
              for i, (path, closed, complete) in enumerate(loops) if closed and complete]
    nearby = defaultdict(list)
    center_tree = cKDTree(xyz[faces].mean(axis=1))
    vertex_tree = cKDTree(xyz)
    hmax = h * settings.max_edge_factor
    records = []
    for ids in groups:
        coord = xyz[ids]
        if np.any((coord <= bmin+1e-8*h) | (coord >= bmax-1e-8*h)) or not curves:
            continue
        near = [np.linalg.norm(v-coord.mean(axis=0), axis=1).min() for _, _, v in curves]
        key = int(np.argmin(near))
        if near[key] <= 2*h:
            nearby[key].extend(ids)
    with tempfile.TemporaryDirectory(dir=out) as temp:
        work = Path(temp)
        replacements, drops = [], set()
        for key, vertices in nearby.items():
            corner, trace, points = curves[key]
            ids = np.unique(vertices)
            center = xyz[ids].mean(axis=0)
            # This operation replaces an entire small corner loop. Preserve large loops.
            if np.linalg.norm(points-center, axis=1).max() > settings.radius_factor*h:
                continue
            distance, ring = vertex_tree.query(points)
            if distance.max() > 1e-8*h or len(np.unique(ring)) != len(ring):
                continue
            adjacent = ((corner+2) % 3+1, corner+1)
            local = np.sort(center_tree.query_ball_point(center, (settings.radius_factor + 4)*h))
            fs, ls = faces[local], labels[local]
            selected = np.isin(fs, np.r_[ring, ids]).any(axis=1)
            best = None
            for layer in range(settings.layers+1):
                if layer:
                    selected = np.isin(fs, np.unique(fs[selected])).any(axis=1)
                selected_ids = local[selected]
                if (not selected.any() or not np.isin(ls[selected], adjacent).all()
                        or any(int(i) in drops for i in selected_ids)):
                    break
                used = np.isin(faces, np.unique(fs[selected])).any(axis=1)
                if not np.isin(np.flatnonzero(used), local).all():
                    break
                try:
                    candidates = list(proposals(xyz, fs, ls, selected, ring, adjacent, hmax))
                except ValueError:
                    continue
                # The outer boundary of the diagnostic neighborhood must survive unchanged.
                local_edges = np.sort(np.vstack([fs[:, [j,(j+1) % 3]] for j in range(3)]), axis=1)
                edges, count = np.unique(local_edges, axis=0, return_counts=True)
                expected_rim = edges[(count == 1) & ~np.isin(edges, ids).all(axis=1)]
                seen = set()
                for t, side, order, longest in candidates:
                    signature = tuple(sorted(tuple(np.roll(tri, -np.argmin(tri)))+(int(label),)
                                             for tri, label in zip(t, side)))
                    if signature in seen:
                        continue
                    seen.add(signature)
                    combined = np.vstack((fs[~selected], t))
                    combined_edges = np.sort(np.vstack([combined[:, [j,(j+1) % 3]] for j in range(3)]), axis=1)
                    edge, number = np.unique(combined_edges, axis=0, return_counts=True)
                    if not np.array_equal(edge[number == 1], expected_rim):
                        continue
                    trial = work/'trial'
                    write_ply(trial.with_suffix('.ply'), xyz, combined, np.r_[ls[~selected], side], corners)
                    result = audit(bin_dir, stages, trial.with_suffix('.ply'), trial, exhaustive=True)
                    if not valid_local(result, len(expected_rim), angle, tolerance):
                        continue
                    score = (result['face_sites_max_projection'], result['projection_p95'], longest)
                    if best is None or score < best[0]:
                        best = (score, selected.copy(), t, side, order, longest, layer, len(expected_rim))
                if best is not None:
                    break
            if best is None:
                continue
            _, selected, t, side, order, longest, layer, expected_edges = best
            discard = local[selected]
            coarse_count = len(t)
            xyz, t, side, fit_levels = fit_replacement(
                bin_dir, stages, work, xyz, fs[~selected], ls[~selected], t, side, corners,
                expected_edges, angle, tolerance, settings.fit_levels)
            drops.update(map(int, discard))
            replacements.append((t, side))
            records.append(dict(center=center.tolist(), corner=int(corner), trace=int(trace), order=order,
                                length_limit=hmax, longest_edge_before_fitting=longest, collar_layers=layer,
                                faces_before_fitting=coarse_count, fit_levels=fit_levels,
                                removed_faces=len(discard), new_faces=len(t), corner_vertices=len(ring)))
        keep = np.ones(len(faces), dtype=bool)
        keep[list(drops)] = False
        ts = np.vstack([faces[keep], *(t for t, _ in replacements)])
        sides = np.concatenate([labels[keep], *(s for _, s in replacements)])
        candidate = work/'candidate.ply'
        write_ply(candidate, xyz, ts, sides, corners)
        final = audit(bin_dir, stages, candidate, out/'final_audit', exhaustive=True)
        old, new = topology(before), topology(final)
        if (not final['intersections_checked'] or final['seam_failures'] or final['seam_chord_failures']
                or final['face_sites_projection_failed'] or final['face_sites_max_projection'] > tolerance
                or any(new[k] > old[k] for k in new if k != 'Euler')):
            raise RuntimeError('The joined candidate failed the full-mesh geometry or topology check')
        candidate.rename(out/'final.ply')
        if replacements:
            write_ply(out/'reconstructed_faces.ply', xyz, np.vstack([t for t,_ in replacements]),
                      np.concatenate([s for _,s in replacements]), (0,0,0))
        complete = all(v == 0 for k, v in new.items() if k != 'Euler')
        report = dict(step=h, repair_max_normal_degrees=angle, max_edge_factor=settings.max_edge_factor,
                      spatial_tolerance=tolerance, pinches_after_rim_pass=len(pinches), repaired=records,
                      topology_before=old, topology_after=new, closed=complete)
        (out/'repair.json').write_text(json.dumps(report, indent=2)+'\n')
        return complete


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path, help='Native closed-range-triangle run with saved stage files')
    parser.add_argument('out', type=Path, help='New results directory')
    parser.add_argument('--bin', type=Path, default=Path(__file__).resolve().parents[3]/'build/src/fiber_surface')
    parser.add_argument('--max-edge-factor', type=float, default=4., help='Maximum new edge length in units of h')
    parser.add_argument('--normal-angle', type=float, help='Maximum normal angle for the replacement neighborhood')
    parser.add_argument('--radius-factor', type=float, default=6., help='Maximum corner-loop radius in units of h')
    parser.add_argument('--fit-levels', type=int, default=2, help='Maximum conforming subdivisions on replacement faces')
    parser.add_argument('--layers', type=int, default=1, help='Additional rings of faces to try around a defect')
    args = parser.parse_args()
    if (not args.max_edge_factor > 0 or not args.radius_factor > 0 or args.layers < 0 or args.fit_levels < 0
            or (args.normal_angle is not None and not 0 < args.normal_angle < 90)):
        parser.error('Expected positive length factors, nonnegative layers, and a normal angle in (0,90)')
    raise SystemExit(0 if reconstruct(args.run, args.out, args.bin, args) else 2)
