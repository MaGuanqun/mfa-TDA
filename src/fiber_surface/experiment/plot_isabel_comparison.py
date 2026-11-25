"""Draw the time-30 Isabel outer-ring cutaway and matched surface details.

Use --volume to select the validated PREFIX_volume.vtk produced by
fiber_volume_mfa for the repaired surface supplied to --ours.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

from join_isabel_refinements import read_ply
from isabel_tetra_section import mesh_section

ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "build"
OUR_MESH = BUILD / "isabel_two_methods/half_fitted/final.ply"
TTK_MESH = BUILD / "fiber-volume-results/isabel_outer_1000x1000x200.3QfVss0i/sides.vtp"
PATCH = np.array([to_rgb(c) for c in ("#70afd5", "#d98ba0", "#b5c66f")])
INK, OPEN, CUT = "#1d303a", "#b02d17", "#f4ca83"
FONT_SIZE = 11.5
X_CUT = 197.
HEIGHT = 2.
SITE = np.array([451.8604, 109.6315, 41.1642])
# Native cell IDs in sides.vtp, checked by the nondegenerate intersection audit.
CROSSING = np.array([1046646, 1867917])
NATIVE_VOLUME = BUILD / "isabel-presentation/native/isabel_volume.vtk"


def read_ttk(path):
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    mesh = reader.GetOutput()
    cells = mesh.GetPolys()
    if not mesh.GetNumberOfPolys() or not np.all(np.diff(vtk_to_numpy(cells.GetOffsetsArray())) == 3):
        raise ValueError(f"Expected a triangle surface in {path}")
    return (vtk_to_numpy(mesh.GetPoints().GetData()).astype(float),
            vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 3).astype(np.int32),
            vtk_to_numpy(mesh.GetCellData().GetArray("EdgeIds")) + 1)


def incidences(mesh):
    """Count on the original mesh IDs, before any cropping or sectioning."""
    points, faces, labels = mesh
    n = len(points)
    keys = []
    for a, b in zip(faces.T, faces[:, [1, 2, 0]].T):
        keys.append(np.minimum(a, b).astype(np.int64)*n + np.maximum(a, b))
    keys = np.concatenate(keys)
    order = keys.argsort()
    unique, start, count = np.unique(keys[order], return_index=True, return_counts=True)
    paired = start[count == 2]
    sides = labels[order % len(faces)]
    ab_bc = ((sides[paired] == 1) & (sides[paired+1] == 2) |
             (sides[paired] == 2) & (sides[paired+1] == 1))
    single = count == 1
    return (np.column_stack(np.divmod(unique[single], n)), sides[start[single]],
            np.column_stack(np.divmod(unique[count == 2][ab_bc], n)))


def unit(v):
    return v / np.linalg.norm(v)


def camera(direction):
    eye = unit(np.asarray(direction, dtype=float))
    right = unit(np.cross([0, 0, 1], eye))
    return np.column_stack((right, np.cross(eye, right), eye))


def section_loops(mesh, cut_x):
    """Cut consistently wound triangles and join segments by their mesh edge IDs."""
    points, faces, _ = mesh
    above = points[:, 0] > cut_x
    faces = faces[above[faces].any(axis=1) & ~above[faces].all(axis=1)]
    start, end = np.empty(len(faces), np.int64), np.empty(len(faces), np.int64)
    positions = {}
    for a, b in zip(faces.T, faces[:, [1, 2, 0]].T):
        crossed = above[a] != above[b]
        up = above[b[crossed]]
        rows = np.flatnonzero(crossed)
        a, b = np.minimum(a[crossed], b[crossed]), np.maximum(a[crossed], b[crossed])
        key = a.astype(np.int64)*len(points) + b
        start[rows[up]], end[rows[~up]] = key[up], key[~up]
        pa, pb = points[a], points[b]
        xyz = pa + ((cut_x-pa[:, 0])/(pb[:, 0]-pa[:, 0]))[:, None]*(pb-pa)
        positions.update(zip(key, xyz))
    successor = dict(zip(start, end))
    if len(successor) != len(faces) or set(successor) != set(end):
        raise ValueError("The section does not consist of closed, oriented contours")
    loops = []
    while successor:
        first = next(iter(successor))
        loop = [first]
        k = successor.pop(first)
        while k != first:
            loop.append(k)
            k = successor.pop(k)
        loops.append(np.array([positions[i] for i in loop + [first]]))
    return loops


def surface_image(mesh, origin, basis, extent, width, *, clip_x=None, wire=False, palette=PATCH):
    """Sample visible triangles with an orthographic depth buffer in grid coordinates.

    Batches bound memory for the 2.6 million TTK triangles. Barycentric
    coordinates also give the distance to the three native edges in the zoom.
    """
    points, faces, labels = mesh
    xmin, xmax, ymin, ymax = extent
    scale = width/(xmax-xmin)
    height = round((ymax-ymin)*scale)
    screen = (points-origin) @ basis
    screen[:, :2] = (screen[:, :2]-[xmin, ymin])*scale - .5
    pixels = np.full((height*width, 3), 255, np.uint8)
    depth = np.full(height*width, -np.inf)
    light = unit(np.array([-.38, .60, .70]))
    for offset in range(0, len(faces), 24000):
        f = faces[offset:offset+24000]
        tri = screen[f]
        lo = np.maximum(np.ceil(tri[:, :, :2].min(axis=1)).astype(int), [0, 0])
        hi = np.minimum(np.floor(tri[:, :, :2].max(axis=1)).astype(int), [width-1, height-1])
        size = np.maximum(hi-lo+1, 0)
        areas = size.prod(axis=1)
        if clip_x is not None:
            areas[points[f, 0].max(axis=1) < clip_x] = 0
        selected = np.flatnonzero(areas)
        # A few triangles fill most of a detail frame. Limit samples per batch.
        groups = np.cumsum(areas[selected]) // 450000
        for ids in np.split(selected, np.flatnonzero(np.diff(groups))+1):
            if not len(ids):
                continue
            t = tri[ids]
            a, b, c = t.transpose(1, 0, 2)
            u, v = b-a, c-a
            normal = np.cross(u, v)
            det = normal[:, 2].copy()
            det[det == 0] = np.nan
            face = np.repeat(np.arange(len(t)), areas[ids])
            start = np.cumsum(areas[ids])-areas[ids]
            k = np.arange(len(face))-start[face]
            x = lo[ids, 0][face] + k % size[ids, 0][face]
            y = lo[ids, 1][face] + k // size[ids, 0][face]
            dx, dy = x-a[face, 0], y-a[face, 1]
            s = (dx*v[face, 1]-dy*v[face, 0])/det[face]
            r = (dy*u[face, 0]-dx*u[face, 1])/det[face]
            z = a[face, 2]+s*u[face, 2]+r*v[face, 2]
            inside = (s >= 0) & (r >= 0) & (s+r <= 1)
            if clip_x is not None:
                world_x = points[f[ids], 0]
                inside &= (world_x[face, 0] + s*(world_x[face, 1]-world_x[face, 0])
                           + r*(world_x[face, 2]-world_x[face, 0])) >= clip_x
            ix = np.flatnonzero(inside)
            location = y[ix]*width+x[ix]
            np.maximum.at(depth, location, z[ix])
            ix = ix[z[ix] >= depth[location]]
            # Undo the pixel scaling before computing two-sided diffuse lighting.
            normal[:, :2] *= scale
            normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1.e-30)
            normal *= np.where(normal[:, 2:] < 0, -1., 1.)
            shade = .69 + .31*np.maximum(0, normal @ light)
            colors = (palette[labels[offset+ids]-1]*shade[:, None])[face[ix]]
            if wire:
                lengths = np.linalg.norm(t[:, [2, 0, 1], :2]-t[:, [1, 2, 0], :2], axis=2)
                altitude = np.abs(det[:, None])/np.maximum(lengths, 1.e-20)
                weights = np.column_stack((1-s[ix]-r[ix], s[ix], r[ix]))
                distance = (weights*altitude[face[ix]]).min(axis=1)
                alpha = .63*np.clip(2.15-distance, 0, 1)
                colors = (1-alpha[:, None])*colors + alpha[:, None]*np.array(to_rgb(INK))
            pixels[y[ix]*width+x[ix]] = np.round(colors*255).astype(np.uint8)
    return np.dstack((pixels.reshape(height, width, 3), 255*np.isfinite(depth.reshape(height, width)).astype(np.uint8)))


def draw(ax, mesh, origin, basis, limits, *, clip_x=None, wire=False, palette=PATCH):
    img = surface_image(mesh, origin, basis, limits, 1350, clip_x=clip_x, wire=wire, palette=palette)
    ax.imshow(img, origin="lower", extent=limits, interpolation="lanczos", aspect="equal", zorder=1)
    ax.set(xlim=limits[:2], ylim=limits[2:], xticks=[], yticks=[])


def edge_faces(faces, edge):
    a, b = edge
    return np.flatnonzero((faces == a).any(axis=1) & (faces == b).any(axis=1))


def transverse_intersection(triangles):
    normals = [unit(np.cross(t[1]-t[0], t[2]-t[0])) for t in triangles]
    direction = unit(np.cross(*normals))
    crossings = []
    for i, triangle in enumerate(triangles):
        dist = (triangle-triangles[1-i, 0]) @ normals[1-i]
        if min(dist.max(), -dist.min()) < 1.e-8:
            raise ValueError("Expected both triangles to penetrate the other's plane")
        hits = list(triangle[np.abs(dist) < 1.e-12])
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if dist[a]*dist[b] < 0:
                hits.append(triangle[a]+dist[a]/(dist[a]-dist[b])*(triangle[b]-triangle[a]))
        crossings.append(np.asarray(hits))
    base = crossings[0][0]
    spans = [(v-base) @ direction for v in crossings]
    ends = np.array([max(s.min() for s in spans), min(s.max() for s in spans)])
    if np.diff(ends).item() <= 0:
        raise ValueError("Selected triangles have no intersection segment")
    return base + ends[:, None]*direction


def closeup_camera(mesh, shared_face_ids, tangent):
    tri = mesh[0][mesh[1][shared_face_ids]]
    normals = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    eye = unit((normals/np.linalg.norm(normals, axis=1)[:, None]).sum(axis=0))
    right = unit(tangent - (tangent @ eye)*eye)
    if np.cross(eye, right)[2] < 0:
        right = -right
    return np.column_stack((right, np.cross(eye, right), eye))


def marked_patch(mesh, center, highlighted):
    points, faces, labels = mesh
    keep = (np.abs(points-center).max(axis=1) < 4)[faces].any(axis=1)
    selected = np.flatnonzero(keep)
    colors = labels[selected].copy()
    for j, face_ids in enumerate(highlighted):
        colors[np.isin(selected, face_ids)] = 4+j
    used, indices = np.unique(faces[keep], return_inverse=True)
    return points[used], indices.reshape(-1, 3), colors


def bounds(xy, pad):
    low, high = xy.min(axis=0)-pad, xy.max(axis=0)+pad
    return low[0], high[0], low[1], high[1]


def axis_indicator(ax, basis, limits):
    # Use the overview camera for all three axes, and anchor the triad in the
    # lower-left margin independently of the chosen data limits.
    xmin, xmax, ymin, ymax = limits
    origin = np.array([xmin+.12*(xmax-xmin), ymin+.055*(ymax-ymin)])
    for label, direction in zip("xyz", basis[:, :2]):
        tip = origin + .073*(xmax-xmin)*direction
        ax.annotate("", xy=tip, xytext=origin, annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", lw=.9, color=INK, mutation_scale=7,
                                    shrinkA=0, shrinkB=0))
        offset = 5*unit(direction)
        ax.annotate(label, tip, xytext=offset, textcoords="offset points", annotation_clip=False,
                    ha="center", va="center", fontsize=FONT_SIZE, color=INK)


def tag(ax, xy, letter, offset):
    ax.scatter(*xy, s=175, facecolor="none", edgecolor=INK, linewidth=1.2, zorder=5)
    ax.annotate(letter, xy, xytext=offset, textcoords="offset points", zorder=6,
                ha="center", va="center", color=INK, fontsize=FONT_SIZE)


def caption(fig, x, y, label):
    fig.text(x, y, label, va="top", ha="center", color=INK, fontsize=FONT_SIZE)


def polygon_area(poly):
    a, b = poly, np.roll(poly, -1, axis=0)
    return .5*abs((a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]).sum())


def section_views(mesh, volume_file):
    loops = section_loops(mesh, X_CUT)
    if len(loops) != 2:
        raise ValueError("The chosen vertical plane must cut two separate limbs of the ring")
    polygons, n_tets = mesh_section(mesh, X_CUT, volume_file)
    # Looking from negative x, increasing y is to the left in the overview.
    loops.sort(key=lambda loop: -loop[:, 1].mean())
    if polygon_area(loops[0][:, 1:]) >= polygon_area(loops[1][:, 1:]):
        raise ValueError("The marked left section must be the smaller one")
    gap = .5*(loops[0][:, 1].min()+loops[1][:, 1].max())
    cuts = [[p for p in polygons if (p[:, 0].mean() > gap) == left] for left in (True, False)]
    for rim, cells in zip(loops, cuts):
        if not np.isclose(sum(map(polygon_area, cells)), polygon_area(rim[:, 1:]), rtol=1e-8):
            raise ValueError("Section cells do not cover the boundary contour")
    return loops, cuts[0], n_tets


def seam_and_intersection(ours, ttk):
    empty, _, shared = incidences(ours)
    if len(empty):
        raise ValueError("Expected the repaired, closed surface")
    unpaired, sides, _ = incidences(ttk)
    selected = np.linalg.norm(ttk[0][unpaired].mean(axis=1)-SITE, axis=1) < .6
    edges = unpaired[selected]
    if sorted(sides[selected]) != [1, 2, 2] or len(np.unique(edges)) != 3:
        raise ValueError("The supplied TTK mesh does not reproduce the open seam")
    seam = ttk[0][np.unique(edges)].mean(axis=0)
    line = transverse_intersection(ttk[0][ttk[1][CROSSING]])
    centers = [seam, line.mean(axis=0)]
    neighbors, cameras = [], []
    for i, center in enumerate(centers):
        midpoints = ours[0][shared].mean(axis=1)
        edge = shared[np.argmin(np.linalg.norm(midpoints-center, axis=1))]
        incident = edge_faces(ours[1], edge)
        if sorted(ours[2][incident]) != [1, 2] or np.linalg.norm(ours[0][edge].mean(axis=0)-center) > .35:
            raise ValueError("Missing the corresponding shared AB/BC edge")
        a, b = ttk[0][edges[sides[selected] == 1][0]] if i == 0 else line
        cameras.append(closeup_camera(ours, incident, b-a))
        neighbors.append(incident)
    return centers, cameras, neighbors, edges, line


def surface_detail(ax, mesh, center, basis, face_ids, segments, mode, wire, half_width):
    # The surrounding mesh is subdued. Strong fills identify the incident faces.
    colors = np.array([to_rgb("#4794b7"), to_rgb("#f7b65e")])
    palette = np.vstack((.075*PATCH + .925, colors))
    groups = [face_ids[mesh[2][face_ids] == label] for label in (1, 2)]
    crop = marked_patch(mesh, center, groups)
    limits = (-half_width, half_width, -.59*half_width, .59*half_width)
    draw(ax, crop, center, basis, limits, wire=wire, palette=palette)
    if mode == "open":
        # Trace the real perimeters of the three singly incident triangles.
        # The three boundary segments nearly coincide at this seam.
        for face in mesh[1][face_ids]:
            perimeter = (mesh[0][face[[0, 1, 2, 0]]]-center) @ basis[:, :2]
            ax.plot(*perimeter.T, color=INK, lw=.55, solid_joinstyle="round", zorder=3)
    if segments is not None:
        for segment in segments:
            xy = (segment-center) @ basis[:, :2]
            ax.plot(*xy.T, color=OPEN, lw=1.2, solid_capstyle="butt", zorder=4)
        if mode == "open":
            nodes = np.unique(segments.reshape(-1, 3), axis=0)
            ax.scatter(*((nodes-center) @ basis[:, :2]).T, color=OPEN, s=6, zorder=5)
    for spine in ax.spines.values():
        spine.set(edgecolor="#a1afb5", linewidth=.55)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--ours", type=Path, default=OUR_MESH)
    cli.add_argument("--ttk", type=Path, default=TTK_MESH)
    cli.add_argument("--volume", type=Path, default=NATIVE_VOLUME, help="Native fiber-volume VTK output")
    cli.add_argument("--output", type=Path, help="PNG/PDF file, or a stem to write both")
    cli.add_argument("--no-wireframe", action="store_true")
    args = cli.parse_args()
    ours, ttk = read_ply(args.ours)[:3], read_ttk(args.ttk)
    centers, cameras, neighbors, edges, line = seam_and_intersection(ours, ttk)
    loops, cells, n_tets = section_views(ours, args.volume)
    left_rim = loops[0]

    plt.rcParams.update({"font.family": "STIXGeneral", "font.size": FONT_SIZE,
                         "font.weight": "normal", "text.color": INK, "pdf.fonttype": 42})
    fig = plt.figure(figsize=(12.1, 4.0), facecolor="white")
    stretch = np.array([1., 1., HEIGHT])
    basis = camera([-.80, -.37, -.48])
    origin = .5*(ours[0].min(axis=0)+ours[0].max(axis=0))*stretch
    def xy(points):
        return (points*stretch-origin) @ basis[:, :2]
    limits = bounds(np.vstack((xy(ours[0]), xy(ttk[0]))), 7)
    top = [(.009+i*.331, .320) for i in range(3)]
    for i, (x, w) in enumerate(top):
        mesh = ttk if i == 0 else ours
        ax = fig.add_axes([x, .470, w, .526])
        ax.set_axis_off()
        draw(ax, (mesh[0]*stretch, *mesh[1:]), origin, basis, limits, clip_x=X_CUT if i == 2 else None)
        if i != 2:
            tag(ax, xy(centers[0]), "df"[i], (7, 12))
            tag(ax, xy(centers[1]), "eg"[i], (3, 13))
        else:
            for loop in loops:
                ax.add_patch(PathPatch(MplPath(xy(loop), closed=True), facecolor=CUT,
                                       edgecolor="#6f441c", linewidth=.65, zorder=2))
            tag(ax, xy(.5*(left_rim.min(axis=0)+left_rim.max(axis=0))), "h", (-13, 5))
        if i == 0:
            axis_indicator(ax, basis, limits)
        caption(fig, x+w/2, .510, ["(a) TTK", "(b) Our method", "(c) Volume cutaway"][i])

    titles = ["(d) TTK · open seam", "(e) TTK · intersecting faces",
              "(f) Ours · closed seam", "(g) Ours · same region", "(h) Left cross section"]
    for i in range(5):
        x, w = .015+.196*i, .185
        ax = fig.add_axes([x, .094, w, .340])
        if i < 4:
            mesh = ttk if i < 2 else ours
            site = i % 2
            center, view = centers[site], cameras[site]
            if i == 0:
                incident = [edge_faces(ttk[1], edge) for edge in edges]
                if any(len(ids) != 1 for ids in incident):
                    raise ValueError("Each highlighted TTK open edge must have one face")
                face_ids, segments, mode = np.concatenate(incident), ttk[0][edges], "open"
            elif i == 1:
                face_ids, segments, mode = CROSSING, line[None, :, :], "crossing"
            else:
                face_ids = neighbors[site]
                segments, mode = None, "shared"
            surface_detail(ax, mesh, center, view, face_ids, segments, mode,
                           not args.no_wireframe, half_width=1.00 if site == 0 else .78)
        else:
            # View normal to the vertical plane, with increasing y to the left.
            scale = np.array([-1., HEIGHT])
            contour = left_rim[:, 1:]*scale
            limits2 = bounds(contour, 1.)
            ax.add_collection(PolyCollection([poly*scale for poly in cells], facecolors="#f9da9d",
                                             edgecolors="#86633a", linewidths=.14))
            ax.plot(*contour.T, color="#553917", lw=.7)
            ax.set(xlim=limits2[:2], ylim=limits2[2:], aspect="equal")
            ax.set_axis_off()
        caption(fig, x+w/2, .078, titles[i])
    name = "isabel-comparison-smooth" if args.no_wireframe else "isabel-comparison"
    output = args.output or ROOT/"paper/figure"/name
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = [output] if output.suffix else [output.with_suffix(ext) for ext in (".png", ".pdf")]
    for path in paths:
        fig.savefig(path, dpi=285 if path.suffix == ".png" else 365)
    plt.close(fig)
    print(f"TTK open edges at {centers[0]}; faces {CROSSING} intersect at {line}. "
          f"Vertical cut crosses two limbs; smaller left section cuts {len(cells)} of {n_tets} tetrahedra.")


if __name__ == "__main__":
    main()
