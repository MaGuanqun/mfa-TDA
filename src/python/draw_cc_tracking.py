import matplotlib as mpl
mpl.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
# Data for the plots
folder = "../../build/src/visualization_result/"

plt.rcParams.update({
    "font.size": 31,
    "mathtext.fontset": "cm",      # Computer Modern
    "font.family": "serif",        # Serif main font
})

fig_size = (8, 5)

def draw_function(ax,xlabel, xlabel_title, num_cc_remove, num_cc, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor,x_label_rotation=0): 
    x_axes = range(len(xlabel))

    # fig, ax = plt.subplots(figsize=fig_size)


    
    ax.plot(num_cc_remove, marker='o', markersize=15, label=r"$\#$CC w/o obstacle")
    ax.plot(num_cc, marker='s', markersize=15, label=r"$\#$CC w/ obstacle")

    ax.set_ylim(50, 1.1*max(max(num_cc_remove), max(num_cc)))
    ax.set_xticks(x_axes)
    # ax.set_xticklabels(xlabel)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)

    ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
    ax.grid(False)

    plt.setp(ax.get_lines(), linewidth=6)
    plt.setp(ax.get_legend().get_lines(), linewidth=6)





def draw_function2(ax,xlabel, xlabel_title, num_cc_remove, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor,x_label_rotation=0): 
    x_axes = range(len(xlabel))
    # fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(num_cc_remove, marker='o', markersize=15, label=r"$\#$CC")

    ax.set_ylim(50, 1.1*max(num_cc_remove))
    ax.set_xticks(x_axes)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)

    ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
    ax.grid(False)


    plt.setp(ax.get_lines(), linewidth=6)
    plt.setp(ax.get_legend().get_lines(), linewidth=6)





def mfa_fig():
    rotation_angle=40
     # Contour data
    step_sizes = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
    xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    # xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    cylinder_cc = [149, 165, 177, 177, 180, 181, 178, 179, 179, 179, 180, 180]
    cylinder_cc_remove_cylinder = [148, 164, 176, 176, 179, 180, 177, 178, 178, 178, 179, 179]
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    vortex_street_cc = [252,257,289,303,306,303,300,298,299,299,299,299]
    vortex_street_cc_remove_cylinder = [251,256,288,302,305,302,299,297,298,298,298,298]
    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,cylinder_cc, 'lower right', (1, 0.1),rotation_angle)
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,vortex_street_cc, 'lower right', (1, 0.1),rotation_angle)
    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'mfa_step_size_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'mfa_step_size_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()

    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'mfa_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'mfa_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()


def inr_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    cylinder_cc = [171,183,211,250,257,276,265,268,256,270,262,275]
    cylinder_cc_remove_cylinder = [167,178,206,239,245,263,251,254,241,256,247,261]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_cc = [118,121,117,124,124,125,124,124,123,123,123,126]
    vortex_street_cc_remove_cylinder = [116,116,112,116,112,118,115,116,113,116,115,115]
    # Draw each plot
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,cylinder_cc, 'lower right', (1, 0.1),rotation_angle)
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,vortex_street_cc, 'lower right', (1, 0.1),rotation_angle)
    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'inr_step_size_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'inr_step_size_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()

    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'inr_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'inr_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()

mfa_fig()
inr_fig()
# s3d_fig()


