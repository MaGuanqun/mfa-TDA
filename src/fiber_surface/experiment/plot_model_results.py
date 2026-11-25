"""Draw the paper figure from the saved full-domain Isabel and upward-wake surfaces."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader

from evaluate_hurricane_isabel import horizontal_cut
from join_isabel_refinements import read_ply

# Side identities follow the analytical teaser. Use darker shades for fine linework.
COLORS = ["#3685b5", "#bd5878", "#839831"]
INK = "#202c34"
REFERENCE = "#a6b0b5"
OPEN = "#b04715"
CENTER = np.array([325.06, 353.24, 42.67])


def top_view(ax, xyz, faces, labels):
    points = xyz[faces]
    normal = np.cross(points[:, 1]-points[:, 0], points[:, 2]-points[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    light = np.array([-.45, .28, .85])
    shade = .58 + .42 * np.abs(normal @ (light / np.linalg.norm(light)))
    rgb = np.asarray([to_rgb(c) for c in COLORS])[labels-1] * shade[:, None]
    split = 292
    for left, right, visible in [(0,split,normal[:,2]>0), (split,499,labels!=1)]:
        ids = np.flatnonzero(visible)
        ids = ids[np.argsort(points[ids, :, 2].mean(axis=1))]
        layer = PolyCollection(points[ids, :, :2], facecolors=rgb[ids],
                               linewidths=0, antialiaseds=False, rasterized=True)
        ax.add_collection(layer)
        layer.set_clip_path(Rectangle((left,0),right-left,499,transform=ax.transData))
    ax.plot([split,split],[323,413],ls=":",lw=.7,color=INK)
    ax.text(split-9,420,"Intact",ha="right",fontsize=7)
    ax.text(split+9,420,"AB hidden",ha="left",fontsize=7)
    ax.add_patch(Circle(CENTER[:2], 8.5, fill=False, edgecolor=INK, linewidth=.9))
    ax.annotate("b", CENTER[:2]+[7,6], xytext=(7,6), textcoords="offset points",
                weight="bold", ha="left", color=INK, arrowprops=dict(arrowstyle="-", lw=.65))
    ax.set(xlim=(0,499), ylim=(0,499), xticks=[0,200,400], yticks=[0,200,400],
           xlabel="$x$ (grid units)", ylabel="$y$ (grid units)", aspect="equal")


def shell_edges(xyz, faces, labels):
    # The isolated shell lies strictly inside this box at every resolution.
    bounds = np.all(np.abs(xyz-CENTER) < [1.5,1.5,.5], axis=1)
    mask = np.all(bounds[faces], axis=1)
    faces, labels = faces[mask], labels[mask]
    edges = np.sort(np.vstack([faces[:, [i,(i+1)%3]] for i in range(3)]), axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    pairs, first, counts = np.unique(edges[order], axis=0, return_index=True, return_counts=True)
    shared = first[counts == 2]
    sides = np.tile(labels, 3)[order]
    seam = (np.minimum(sides[shared],sides[shared+1])==1) & (np.maximum(sides[shared],sides[shared+1])==2)
    return xyz[pairs[counts==2][seam], :2], xyz[pairs[counts==1], :2]


def zoom(ax, mesh, ref, title):
    seam, exposed = shell_edges(*mesh)
    for segs, color, width, style in [(ref,REFERENCE,2.6,"solid"), (seam,INK,1.3,"solid"),
                                      (exposed,OPEN,1.2,"dashed")]:
        ax.add_collection(LineCollection(segs, colors=color, linewidths=width,
                                         linestyles=style, capstyle="round"))
    ax.set(xlim=(324.15,325.43), ylim=(352.74,353.86), aspect="equal")
    ax.set_axis_off()
    ax.set_title(title, fontsize=8, y=.95, pad=1)
    ax.plot([324.23,324.48],[352.77,352.77], color=INK, lw=1)
    ax.text(324.355,352.785,"0.25", ha="center", va="bottom", fontsize=6)


def wake(ax, file, title, show_x):
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(str(file))
    reader.Update()
    segments, labels = horizontal_cut(reader.GetOutput(), 14.1)
    if set(np.unique(labels)) != {1,2,3}:
        raise ValueError(f"Missing range side in the section of {file}")
    ax.add_collection(LineCollection(segments, colors=np.asarray(COLORS)[labels-1],
                                     linewidths=1.0, capstyle="round"))
    ax.set(xlim=(-.5,7.5), ylim=(-.5,.5), yticks=[-.5,.5], xticks=[0,2,4,6],
           ylabel="$y$", aspect="equal")
    ax.tick_params(axis="x", labelbottom=show_x)
    if show_x:
        ax.set_xlabel("$x$", labelpad=-7, loc="right")
    ax.text(.985,.90,title, transform=ax.transAxes, ha="right", va="top", fontsize=7)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--build", type=Path, default=Path(__file__).resolve().parents[3]/"build")
    cli.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[3]/"paper/figure/model-results")
    args = cli.parse_args()
    half = args.build/"isabel_half_study"
    studies = args.build/"isabel_two_methods"
    wake_dir = args.build/"fiber-volume-results/vortex_street_up"
    fitted = read_ply(studies/"half_fitted/final.ply")[:3]
    raw_paths = list(half.glob("hurricane_isabel_t30_outer_h0.5*/volume_boundary.ply"))
    if len(raw_paths)!=1:
        raise ValueError(f"Expected one native half-step boundary in {half}")
    raw = read_ply(raw_paths[0])[:3]
    ref, _ = shell_edges(*read_ply(half/"site325/ref.03125.ply")[:3])

    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":7.5, "pdf.fonttype":42,
                         "axes.linewidth":.45, "xtick.major.width":.45,"ytick.major.width":.45,
                         "xtick.major.size":2.0,"ytick.major.size":2.0,"axes.labelpad":1.5})
    fig = plt.figure(figsize=(7.05,3.28))
    def heading(x,y,text):
        fig.text(x,y,text,ha="left",va="top",weight="bold",fontsize=9,color=INK)
    heading(.018,.988,"(a)  Hurricane Isabel · repaired surface")
    heading(.482,.988,"(b)  Isolated shell · corner B")
    top_view(fig.add_axes([.074,.225,.368,.695]), *fitted)
    zoom(fig.add_axes([.520,.567,.190,.344]), raw, ref, "Before repair")
    zoom(fig.add_axes([.756,.567,.190,.344]), fitted, ref, "After repair")
    zoom_keys = [Line2D([],[],color=c,lw=w,ls=s,label=t) for c,w,s,t in
                 [(REFERENCE,2.6,"-","Reference"),(INK,1.3,"-","Shared seam"),(OPEN,1.2,"--","Open edge")]]
    fig.legend(handles=zoom_keys,loc="center",bbox_to_anchor=(.737,.535),
               ncol=3,frameon=False,fontsize=6.8,handlelength=1.5,columnspacing=1.)
    heading(.482,.475,"(c)  Vortex street · upward flow at $t=14.1$")
    wake(fig.add_axes([.531,.250,.450,.138]), wake_dir/"h005/side_patches.vtp", "Our mesh", False)
    wake(fig.add_axes([.531,.082,.450,.138]), wake_dir/"ttk/grid0025/sides.vtp", "TTK", True)
    fig.legend(handles=[Patch(facecolor=c,label=t) for c,t in zip(COLORS,("AB","BC","CA"))],
               loc="center", bbox_to_anchor=(.252,.078), ncol=3, frameon=False,
               handlelength=1.,columnspacing=1.5, title="Range edges",title_fontsize=7.2)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","png"):
        fig.savefig(args.output.with_suffix("."+ext),dpi=400 if ext=="pdf" else 220)
    plt.close(fig)


if __name__ == "__main__":
    main()
