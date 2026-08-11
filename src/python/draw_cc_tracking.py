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
                  legend_loc, bbox_to_anchor,x_label_rotation=0, ylim=50, label1=r"$\#$CC w/ obstacle", label2=r"$\#$CC w/o obstacle", y_max=None): 
    x_axes = range(len(xlabel))

    # fig, ax = plt.subplots(figsize=fig_size)


    
    ax.plot(num_cc, marker='s', markersize=20, label=label1)
    ax.plot(num_cc_remove, marker='o', markersize=15, label=label2)

    if y_max is not None:
        ax.set_ylim(ylim, y_max)
    else:
        ax.set_ylim(ylim, 1.1*max(max(num_cc_remove), max(num_cc)))
    ax.set_xticks(x_axes)
    # ax.set_xticklabels(xlabel)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)

    ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
    ax.grid(False)

    plt.setp(ax.get_lines(), linewidth=6)
    plt.setp(ax.get_legend().get_lines(), linewidth=6)


def draw_function_3_curves(ax,xlabel, xlabel_title, num_cc_remove, num_cc, third_curve, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor,x_label_rotation=0, ylim=50, label1=r"$\#$CC w/ obstacle", label2=r"$\#$CC w/o obstacle", label3="third curve", y_max=None): 
    x_axes = range(len(xlabel))

    # fig, ax = plt.subplots(figsize=fig_size)

    lw_first = 9
    ax.plot(num_cc, marker='s', markersize=20, label=label1, linewidth=lw_first, color="C0")
    ax.plot(num_cc_remove, marker='o', markersize=15, label=label2, color="green")
    ax.plot(third_curve, marker='^', markersize=15, label=label3, color="C1")
    
    if y_max is not None:
        ax.set_ylim(ylim, y_max)
    else:
        ax.set_ylim(ylim, 1.1*max(max(num_cc_remove), max(num_cc)))
    ax.set_xticks(x_axes)
    # ax.set_xticklabels(xlabel)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)

    ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
    ax.grid(False)



def draw_function2(ax,xlabel, xlabel_title, num_cc_remove, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor,x_label_rotation=0, ylim=50,label1=r"$\#$CC",y_max=None): 
    x_axes = range(len(xlabel))
    # fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(num_cc_remove, marker='o', markersize=15, label=label1)

    if y_max is not None:
        ax.set_ylim(ylim, y_max)
    else:
        ax.set_ylim(ylim, 1.1*max(num_cc_remove))
    ax.set_xticks(x_axes)
    
    ax.set_xticklabels(xlabel, rotation=x_label_rotation)
    
    ax.set_xlabel(xlabel_title)
    plt.setp(ax.get_lines(), linewidth=6)
    if label1 is not None:
        ax.legend(loc=legend_loc, bbox_to_anchor=bbox_to_anchor)
        plt.setp(ax.get_legend().get_lines(), linewidth=6)
        
    ax.grid(False)






def mfa_fig():
    rotation_angle=40
     # Contour data
    step_sizes = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
    xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    # xlabel_contour = [r"$l$/2", r"$l$/4", r"$l$/8", r"$l$/16", r"$l$/32", r"$l$/64", r"$l$/128", r"$l$/256"]    
    cylinder_cc = [238, 279, 291, 304, 308, 301, 302, 302, 302, 302, 302, 302]
    cylinder_cc_remove_cylinder = [237, 278, 290, 303, 307, 300, 301, 301, 301, 301, 301, 301]
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    vortex_street_cc = [245, 261, 292, 299, 302, 301, 300, 298, 299, 300, 300, 300]
    vortex_street_cc_remove_cylinder = [244, 260, 291, 298, 301, 300, 299, 297, 298, 299, 299, 299]
    # Draw each plot
    
    
    vortex_street_loop=[5, 5, 11, 18, 22, 36, 36, 40, 42, 50, 55, 53]
    vortex_street_loop_remove_cylinder=[5, 5, 11, 18, 22, 36, 36, 40, 42, 50, 55, 53]
    
    cylinder_loop=[2, 16, 23, 39, 41, 53, 59, 65, 72, 72, 80, 81]
    cylinder_loop_remove_cylinder=[2, 16, 23, 39, 41, 53, 59, 65, 72, 72, 80, 81]


    vortex_cc = [7702, 10551, 11562, 10891, 10239, 9573, 9351, 9149, 9060, 8975, 8953, 8936]
    vortex_loop = [59, 203, 557, 1764, 2906, 4247, 5047, 5668, 5834, 6068, 6763, 6819]

    
    label1=r"$\#$Loop w/ obstacle"
    label2=r"$\#$Loop w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,cylinder_cc, 'lower right', (1.0, 1),rotation_angle,150)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, cylinder_loop_remove_cylinder, cylinder_loop, 'lower right', (1, 1.0),rotation_angle,0, label1, label2)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_mfa_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_mfa_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]}) 
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,vortex_street_cc, 'lower right', (1, 1.0),rotation_angle,150)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_loop_remove_cylinder,vortex_street_loop, 'lower right', (1, 1.0),rotation_angle,0, label1, label2)
    
    # Adjust layout and save the figure
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_cylinder_mfa_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_cylinder_mfa_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    label5=r"$\#$Loop"
    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]}) 
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_cc,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_loop,'lower right', (1, 0.1), rotation_angle,0,label5)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_mfa_curve.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_mfa_curve.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    rotation_angle=60
    
    fig, axes = plt.subplots(1, 3, figsize=(20, 4),gridspec_kw={'width_ratios': [5.5, 5.5, 5.5]})  # 3 subplots in a row
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,150)
    
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,150)
    
    draw_function2(axes[2], xlabel_contour, xlabel_title_contour, vortex_cc,'lower right', (1, 0.1), rotation_angle,0)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'mfa_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'mfa_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    


def inr_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    cylinder_cc = [203, 254, 314, 305, 312, 333, 348, 328, 334, 333, 339, 328]
    cylinder_cc_remove_cylinder = [189, 229, 263, 288, 285, 290, 284, 272, 271, 278, 277, 281]
    # cylinder_cc = [202, 268, 302, 326, 350, 378, 385, 390, 398, 407, 399, 408]
    # cylinder_cc_remove_cylinder = [191, 240, 274, 302, 308, 313, 301, 294, 297, 307, 303, 318]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    vortex_street_cc = [133, 127, 137, 134, 145, 141, 151, 149, 149, 149, 155, 152]
    vortex_street_cc_remove_cylinder = [125, 123, 122, 123, 122, 123, 128, 125, 125, 125, 125, 122]
    
    
    # vortex_street_cc = [128, 135, 141, 149, 157, 162, 165, 163, 165, 167, 166, 173]
    # vortex_street_cc_remove_cylinder = [125, 124, 124, 124, 125, 127, 129, 128, 125, 127, 125, 126]
    
    
    vortex_street_loop = [3, 4, 8, 6, 7, 8, 6, 3, 3, 4, 3, 4]
    vortex_street_loop_remove_cylinder=[1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 3]
    
    # vortex_street_loop=[3, 5, 7, 4, 2, 2, 1, 1, 2, 1, 1, 1]
    # vortex_street_loop_remove_cylinder=[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    
    cylinder_loop = [10, 19, 23, 25, 32, 24, 27, 36, 48, 50, 53, 53]
    
    cylinder_loop_remove_cylinder = [9, 17, 18, 13, 17, 23, 26, 35, 47, 49, 52, 51]
    
    # cylinder_loop=[7, 17, 19, 18, 23, 24, 20, 29, 39, 42, 41, 41]
    # cylinder_loop_remove_cylinder=[7, 16, 17, 14, 16, 17, 20, 29, 39, 42, 41, 41]
    
    label1=r"$\#$Loop w/ obstacle"
    label2=r"$\#$Loop w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,cylinder_cc, 'lower right', (1, 1.0),rotation_angle,0)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, cylinder_loop_remove_cylinder, cylinder_loop, 'lower right', (1.0, 1.0),rotation_angle,0, label1, label2)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_inr_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_inr_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
        
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,vortex_street_cc, 'lower right', (1, 1.0),rotation_angle,0)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_loop_remove_cylinder,vortex_street_loop, 'lower right',(1,1.0),rotation_angle,0, label1, label2)
    # Adjust layout and save the figure

    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_cylinder_inr_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_cylinder_inr_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    

    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,0)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'inr_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'inr_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    





def inr_matching_score_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    
    cylinder_discrete_score_with_obstacle = [0.909, 0.843, 0.824, 0.762, 0.727, 0.690, 0.619, 0.579, 0.519, 0.477, 0.442, 0.400]
    cylinder_discrete_score_without_obstacle =[0.895, 0.837, 0.821, 0.771, 0.736, 0.703, 0.669, 0.655, 0.639, 0.626, 0.621, 0.607]
    
    cylinder_implicit_score_with_obstacle=[0.740, 0.779, 0.803, 0.843, 0.878, 0.899, 0.915, 0.933, 0.941, 0.957, 0.961, 0.965]
    
    cylinder_implicit_score_without_obstacle=[0.745, 0.780, 0.825, 0.867, 0.899, 0.928, 0.942, 0.951, 0.956, 0.963, 0.965, 0.969]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_discrete_score_with_obstacle = [0.930, 0.934, 0.918, 0.866, 0.823, 0.736, 0.660, 0.583, 0.528, 0.478, 0.433, 0.380]
    vortex_street_discrete_score_without_obstacle =  [0.932, 0.939, 0.948, 0.935, 0.922, 0.898, 0.879, 0.849, 0.829, 0.797, 0.759, 0.717]
    
    
    vortex_street_implicit_score_with_obstacle=[0.923, 0.924, 0.946, 0.956, 0.960, 0.973, 0.979, 0.984, 0.985, 0.986, 0.988, 0.988]
        
        
    vortex_street_implicit_score_without_obstacle = [0.945, 0.957, 0.964, 0.969, 0.975, 0.980, 0.984, 0.987, 0.987, 0.989, 0.990, 0.989]
    
    label1=r"w/ obstacle"
    label2=r"w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_discrete_score_without_obstacle,cylinder_discrete_score_with_obstacle, 'lower right', (0.62, 0.0),rotation_angle,0,label1,label2, 1.1)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, cylinder_implicit_score_without_obstacle, cylinder_implicit_score_with_obstacle, 'lower right', (1.0, 0.0),rotation_angle,0, label1, label2,1.1)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_inr_score.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_inr_score.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
        
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_discrete_score_without_obstacle,vortex_street_discrete_score_with_obstacle, 'lower right', (0.62, 0.0),rotation_angle,0, label1, label2,1.1)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_implicit_score_without_obstacle,vortex_street_implicit_score_with_obstacle, 'lower right',(1,0.0),rotation_angle,0, label1, label2,1.1)
    # Adjust layout and save the figure

    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_cylinder_inr_score.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_cylinder_inr_score.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    

def mfa_matching_score_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
    xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    
    
    cylinder_discrete_score_with_obstacle = [0.941, 0.952, 0.958, 0.971, 0.975, 0.984, 0.987, 0.989, 0.988, 0.989, 0.990, 0.992]
    cylinder_discrete_score_without_obstacle =[0.936, 0.948, 0.953, 0.968, 0.973, 0.981, 0.985, 0.988, 0.987, 0.988, 0.989, 0.991]
    
    cylinder_implicit_score_with_obstacle= [0.864, 0.918, 0.922, 0.948, 0.962, 0.979, 0.983, 0.988, 0.989, 0.990, 0.990, 0.991]
    
    cylinder_implicit_score_without_obstacle=[0.858, 0.914, 0.919, 0.948, 0.960, 0.977, 0.982, 0.988, 0.987, 0.989, 0.990, 0.990]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_discrete_score_with_obstacle =[0.950, 0.909, 0.950, 0.947, 0.969, 0.976, 0.974, 0.977, 0.978, 0.982, 0.980, 0.980]
    vortex_street_discrete_score_without_obstacle =  [0.950, 0.908, 0.950, 0.947, 0.969, 0.976, 0.974, 0.977, 0.978, 0.982, 0.980, 0.980]
    
    
    vortex_street_implicit_score_with_obstacle=[0.709, 0.852, 0.895, 0.939, 0.951, 0.968, 0.974, 0.981, 0.983, 0.888, 0.987, 0.988]
        
        
    vortex_street_implicit_score_without_obstacle =[0.707, 0.851, 0.895, 0.939, 0.951, 0.967, 0.974, 0.980, 0.983, 0.888, 0.987, 0.988]
    
    label1=r"w/ obstacle"
    label2=r"w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_discrete_score_without_obstacle,cylinder_discrete_score_with_obstacle, 'lower right', (1.0, 0.0),rotation_angle,0,label1,label2, 1.1)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, cylinder_implicit_score_without_obstacle, cylinder_implicit_score_with_obstacle, 'lower right', (1.0, 0.0),rotation_angle,0, label1, label2,1.1)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_mfa_score.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_mfa_score.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
        
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_discrete_score_without_obstacle,vortex_street_discrete_score_with_obstacle, 'lower right', (1, 0.0),rotation_angle,0, label1, label2,1.1)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_implicit_score_without_obstacle,vortex_street_implicit_score_with_obstacle, 'lower right',(1,0.0),rotation_angle,0, label1, label2,1.1)
    # Adjust layout and save the figure

    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_cylinder_mfa_score.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_cylinder_mfa_score.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    rotation_angle=0
     # Contour data
    step_sizes = [2, 3, 4, 6, 8, 12, 16]
    xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    
    
    vortex_mfa_discrete_score= [0.773, 0.682, 0.700, 0.769, 0.796, 0.809, 0.812]
    vortex_mfa_continuous_score=[0.692, 0.673, 0.674, 0.684, 0.690, 0.694, 0.695]
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_mfa_discrete_score,'lower right', (1, 0.1), rotation_angle,0, None,1)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_mfa_continuous_score,'lower right', (1, 0.1), rotation_angle,0,None,1)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_mfa_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_mfa_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
def gradient_norm_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    
    vortex_inr_mean = [1.522e-11, 1.573e-11, 1.477e-11, 1.449e-11, 1.382e-11, 1.211e-11, 1.085e-11, 9.722e-12, 8.501e-12, 7.478e-12, 6.927e-12, 6.303e-12]

    vortex_inr_max = [9.976e-11, 9.990e-11, 9.993e-11, 9.996e-11, 9.989e-11, 9.999e-11, 9.998e-11, 9.980e-11, 1.000e-10, 9.997e-11, 9.998e-11, 1.000e-10]

    cylinder_inr_mean =[1.130e-11, 1.171e-11, 1.157e-11, 1.292e-11, 1.229e-11, 1.288e-11, 1.289e-11, 1.299e-11, 1.317e-11, 1.320e-11, 1.310e-11, 1.329e-11]
    
    cylinder_inr_max = [9.989e-11, 9.997e-11, 9.999e-11, 9.997e-11, 9.989e-11, 9.999e-11, 9.998e-11, 9.999e-11, 1.000e-10, 9.999e-11, 9.999e-11, 1.000e-10]

    vortex_mfa_mean = [6.385e-12, 5.984e-12, 6.426e-12, 6.259e-12, 4.637e-12, 1.693e-12, 2.158e-12, 6.181e-12, 1.467e-11, 2.540e-11, 2.808e-11, 2.708e-11]

    vortex_mfa_max = [9.966e-11, 9.998e-11, 9.978e-11, 9.986e-11, 9.991e-11, 9.984e-11, 9.984e-11, 9.997e-11, 9.999e-11, 1.000e-10, 1.000e-10, 1.000e-10]

    cylinder_mfa_mean = [4.164e-12, 5.380e-12, 9.998e-12, 1.643e-11, 1.751e-11, 1.402e-11, 1.304e-11, 1.401e-11, 1.646e-11, 1.929e-11, 2.029e-11, 2.012e-11]

    cylinder_mfa_max =[9.993e-11, 9.991e-11, 9.994e-11, 9.998e-11, 1.000e-10, 9.992e-11, 9.997e-11, 9.998e-11, 9.999e-11, 1.000e-10, 1.000e-10, 1.000e-10]
    
    
    vortex_4d_mfa_mean=[6.264e-12, 6.004e-12, 5.803e-12, 7.980e-12, 8.869e-12, 2.945e-11, 3.695e-12, 2.281e-12, 4.822e-12, 1.797e-11, 6.141e-12, 3.577e-11]
    
    vortex_4d_mfa_max=[9.998e-11, 1.000e-10, 9.998e-11, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10, 1.000e-10]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 1, figsize=(8.5, 5))  # 3 subplots in a row

    
    label1=r"Max"
    label2=r"Mean"
    
    
    draw_function(axes, xlabel_contour, xlabel_title_contour, vortex_inr_mean,vortex_inr_max, 'lower right', (0.6, 0.2),rotation_angle,0,label1,label2, 1.1e-10)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_inr_grad_norm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_inr_grad_norm.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    
    fig, axes = plt.subplots(1, 1, figsize=(8.5, 5))  # 3 subplots in a row

    draw_function(axes, xlabel_contour, xlabel_title_contour, cylinder_inr_mean,cylinder_inr_max, 'lower right', (0.6, 0.2),rotation_angle,0,label1,label2, 1.1e-10)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_inr_grad_norm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_inr_grad_norm.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    step_sizes = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
    xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    
    
    
    
        
    fig, axes = plt.subplots(1, 1, figsize=(8.5, 5))  # 3 subplots in a row
    
    draw_function(axes, xlabel_contour, xlabel_title_contour, vortex_mfa_mean,vortex_mfa_max, 'lower right', (0.6, 0.2),rotation_angle,0,label1,label2, 1.1e-10)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_mfa_grad_norm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_mfa_grad_norm.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 1, figsize=(8.5, 5))  # 3 subplots in a row

    draw_function(axes, xlabel_contour, xlabel_title_contour, cylinder_mfa_mean,cylinder_mfa_max, 'lower right', (0.6, 0.2),rotation_angle,0,label1,label2, 1.1e-10)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_mfa_grad_norm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_mfa_grad_norm.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    fig, axes = plt.subplots(1, 1, figsize=(8.5, 5))  # 3 subplots in a row
    draw_function(axes, xlabel_contour, xlabel_title_contour, vortex_4d_mfa_mean,vortex_4d_mfa_max, 'lower right', (0.6, 0.2),rotation_angle,0,label1,label2, 1.1e-10)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_mfa_grad_norm.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_mfa_grad_norm.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
def closed_form_fig():
    
    rotation_angle=70
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    rotation_cc = [9,9,9,9,9,9,9,9,9,9,9,9]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    
    potential_cc = [1,1,1,1,1,1,1,1,1,1,1,1]
       # Draw each plot
    fig, axes = plt.subplots(1, 3, figsize=(17, 4),gridspec_kw={'width_ratios': [5.5, 5.5, 5.5]})  # 3 subplots in a row
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, potential_cc,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, rotation_cc,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[2], xlabel_contour, xlabel_title_contour, potential_cc,'lower right', (1, 0.1), rotation_angle,0)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'closed_form_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'closd_form_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()


# def draw_4d_fig():
#     rotation_angle=40
#     step_sizes = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
#     xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
#     vortex_cc = [7702, 10551, 11562, 10891, 10239, 9573, 9351, 9149, 9060, 8975, 8953, 8934]
#     vortex_loop = [59, 203, 557, 1764, 2906, 4247, 5047, 5668, 5834, 6068, 6763, 6212]
#     xlabel_title_contour = r"Step size $s$"
#     fig, axes = plt.subplots(1, 1, figsize=(8, 4))  # 3 subplots in a row
#     draw_function2(axes, xlabel_contour, xlabel_title_contour, vortex_cc, 'lower right', (1, 0.1),rotation_angle,0)
#     fig.tight_layout(pad=0.1)
#     fig.savefig(folder + 'vortex_4d.png', bbox_inches='tight', pad_inches=0)
#     fig.savefig(folder + 'vortex_4d.pdf', bbox_inches='tight', pad_inches=0)
#     fig.clear()
#     plt.close(fig)
    
def draw_start_points_with_seeds():
    rotation_angle=40
    xlabel_contour = [20, 30, 40, 60, 80, 120, 160, 240, 320]
    # xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    boussinesq_inr = [111,144,164,172,175,175,175,175,175]
    vortex_street_inr=[84, 87,87,89,90,92,92,91,91]
    xlabel_title_contour = "Number of Seeds"
    
    quartic_potential_2d =[1,1,1,1,1,1,1,1,1]
    quartic_potential_3d = [1,1,1,1,1,1,1,1,1]
    rotating_quartic_multiwell =[9,9,9,9,9,9,9,9,9]
    
    label1=r"Quartic Potential 2D"
    label2=r"Quartic Rotation"
    label3=r"Quartic Potential 3D"
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 7),gridspec_kw={'width_ratios': [5.5, 5.5]}) 
    
    draw_function_3_curves(axes[0], xlabel_contour, xlabel_title_contour,rotating_quartic_multiwell, quartic_potential_2d,quartic_potential_3d, 'lower right',(1,1.0),rotation_angle,0, label1, label2,label3)
        
    label1=r"Vortex Street INR"
    label2=r"Heated Cylinder INR"
    
        
    draw_function(axes[1], xlabel_contour, xlabel_title_contour,boussinesq_inr, vortex_street_inr,'lower right',(1,1.0),rotation_angle,0, label1, label2)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'start_points_number_with_seeds.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'start_points_number_with_seeds.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    plt.close(fig)
    
    
        
def draw_start_points_with_seeds_analytic():
    rotation_angle=40
    xlabel_contour = [20, 30, 40, 60, 80, 120, 160, 240, 320]
    # xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    quartic_potential_2d =[1,1,1,1,1,1,1,1,1]
    quartic_potential_3d = [1,1,1,1,1,1,1,1,1]
    rotating_quartic_multiwell =[9,9,9,9,9,9,9,9,9]

    xlabel_title_contour = "Number of Seeds"
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))  # 3 subplots in a row
    label1=r"Quartic Potential 2D"
    label2=r"Quartic Rotation"
    label3=r"Quartic Potential 3D"
    
    draw_function_3_curves(axes, xlabel_contour, xlabel_title_contour,rotating_quartic_multiwell, quartic_potential_2d,quartic_potential_3d, 'lower right',(1,1.0),rotation_angle,0, label1, label2,label3)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'analytic_start_points_number_with_seeds.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'analytic_start_points_number_with_seeds.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    plt.close(fig)
    


def draw_start_points_with_epsilon():
    rotation_angle=40
    xlabel_contour = ["$e^{-6}$", "$e^{-7}$", "$e^{-8}$", "$e^{-9}$", "$e^{-10}$","$e^{-11}$","$e^{-12}$", "$e^{-13}$"]
    # xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    boussinesq_inr = [147,147,147,147,147,147,147,147]
    vortex_street_inr=[83,83,83,83,83,83,83,83]
    
    boussinesq_mfa = [206,206,206,206,206,206,206,206]
    vortex_street_mfa=[158,158,158,158,158,158,158,158]
    
    vortex_mfa=[7082,7079,7079,7079,7079,7079,7079,7079]
    xlabel_title_contour = r"$\epsilon$"
    fig, axes = plt.subplots(1, 3, figsize=(22, 6),gridspec_kw={'width_ratios': [6.5, 6.5, 6.5]})  # 3 subplots in a row
    
    lable1=r"Vortex Street MFA"
    lable2=r"Heated Cylinder MFA"
    draw_function(axes[0], xlabel_contour, xlabel_title_contour,boussinesq_mfa,vortex_street_mfa, 'lower right',(1,1.0),rotation_angle,50, lable1, lable2, 250)
    
    lable1=r"Vortex MFA"
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_mfa,'lower right', (1, 1.0), rotation_angle,5000,lable1,8000)
    
    label1=r"Vortex Street INR"
    label2=r"Heated Cylinder INR"
    draw_function(axes[2], xlabel_contour, xlabel_title_contour,boussinesq_inr,vortex_street_inr, 'lower right',(1,1.0),rotation_angle,30, label1, label2,200)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'start_points_number_with_epsilon.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'start_points_number_with_epsilon.pdf', bbox_inches='tight', pad_inches=0)



def draw_degenerate_points_with_epsilon_k():
    rotation_angle=40
    xlabel_contour = ["$e^{-2}$", "$e^{-3}$", "$e^{-4}$", "$e^{-5}$", "$e^{-6}$","$e^{-7}$","$e^{-8}$"]
    # xlabel_contour = [rf"$l/{s}$" for s in step_sizes]
    boussinesq_mfa = [265,265,265,265,265,265,265]
    vortex_street_mfa=[221,221,221,221,221,221,221]
    
    boussinesq_inr = [693,693,693,693,693,693,693]
    vortex_street_inr=[195,195,195,195,195,195,195]
    
    vortex_mfa=[30398,30398,30398,30398,30398,30398,30398]
    xlabel_title_contour = r"$\epsilon_k$"
    fig, axes = plt.subplots(1, 3, figsize=(22, 6),gridspec_kw={'width_ratios': [6.5, 6.5, 6.5]})  # 3 subplots in a row
    
    lable1=r"Vortex Street MFA"
    lable2=r"Heated Cylinder MFA"
    draw_function(axes[0], xlabel_contour, xlabel_title_contour,boussinesq_mfa,vortex_street_mfa, 'lower right',(1,1.0),rotation_angle,100, lable1, lable2, 300)
    
    lable1=r"Vortex MFA"
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, vortex_mfa,'lower right', (1, 1.0), rotation_angle,20000,lable1,40000)
    
    label1=r"Vortex Street INR"
    label2=r"Heated Cylinder INR"
    draw_function(axes[2], xlabel_contour, xlabel_title_contour,boussinesq_inr,vortex_street_inr, 'lower right',(1,1.0),rotation_angle,50, label1, label2,750)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'degenerate_points_number_with_epsilon_k.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'degenerate_points_number_with_epsilon_k.pdf', bbox_inches='tight', pad_inches=0)
    
# mfa_fig() 
# inr_fig()
# closed_form_fig()
# inr_matching_score_fig() 

# mfa_matching_score_fig()
# draw_4d_fig()
# gradient_norm_fig()
# draw_start_points_with_seeds()
draw_start_points_with_epsilon()
# draw_start_points_with_seeds_analytic()
# draw_degenerate_points_with_epsilon_k()