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





def draw_function2(ax,xlabel, xlabel_title, num_cc_remove, # fig_name, pdf_name, 
                  legend_loc, bbox_to_anchor,x_label_rotation=0, ylim=50): 
    x_axes = range(len(xlabel))
    # fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(num_cc_remove, marker='o', markersize=15, label=r"$\#$CC")

    ax.set_ylim(ylim, 1.1*max(num_cc_remove))
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
    cylinder_cc = [242,276,296,303,307,300,301,301,301,301,301,302]
    cylinder_cc_remove_cylinder = [241,275,295,302,306,299,300,300,300,300,300,301]
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    vortex_street_cc = [252,257,289,303,306,303,300,298,299,299,299,299]
    vortex_street_cc_remove_cylinder = [251,256,288,302,305,302,299,297,298,298,298,298]
    # Draw each plot
    
    
    vortex_street_loop=[12, 10, 18, 22, 31, 42, 55, 51, 54, 66, 70, 74]
    vortex_street_loop_remove_cylinder=[12, 10, 18, 22, 31, 42, 55, 51, 54, 66, 70, 74]
    
    cylinder_loop=[2, 11, 19, 26, 33, 39, 57, 53, 51, 70, 60, 59]
    cylinder_loop_remove_cylinder=[2, 11, 19, 26, 33, 39, 57, 53, 51, 70, 60, 59]


    
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
    
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,150)
    
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle,150)
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
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_cc = [118,121,117,124,124,125,124,124,123,123,123,126]
    vortex_street_cc_remove_cylinder = [116,116,112,116,112,118,115,116,113,116,115,115]
    
        
    vortex_street_loop=[4, 4, 3, 3, 7, 2, 3, 4, 6, 4, 6, 2]
    vortex_street_loop_remove_cylinder=[3, 3, 1, 3, 5, 2, 2, 2, 5, 4, 5, 2]
    
    cylinder_loop=[26, 29, 29, 19, 17, 13, 19, 18, 23, 21, 28, 31]
    cylinder_loop_remove_cylinder=[24, 27, 26, 19, 17, 13, 19, 18, 23, 21, 28, 31]
    
    label1=r"$\#$Loop w/ obstacle"
    label2=r"$\#$Loop w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,cylinder_cc, 'lower right', (1, 1.0),rotation_angle)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, cylinder_loop_remove_cylinder, cylinder_loop, 'lower right', (1.0, 1.0),rotation_angle,0, label1, label2)
    
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'heated_cylinder_inr_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'heated_cylinder_inr_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(17, 6),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
        
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,vortex_street_cc, 'lower right', (1, 1.0),rotation_angle)
    
    draw_function(axes[1], xlabel_contour, xlabel_title_contour, vortex_street_loop_remove_cylinder,vortex_street_loop, 'lower right',(1,1.0),rotation_angle,0, label1, label2)
    # Adjust layout and save the figure

    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'vortex_street_cylinder_inr_two_curves.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'vortex_street_cylinder_inr_two_curves.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    

    fig, axes = plt.subplots(1, 2, figsize=(17, 4),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, cylinder_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, vortex_street_cc_remove_cylinder,'lower right', (1, 0.1), rotation_angle)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'inr_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'inr_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()
    





def inr_matching_score_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    
    cylinder_discrete_score_with_obstacle = [0.831, 0.671, 0.617, 0.525, 0.483, 0.427, 0.388, 0.348, 0.311, 0.284, 0.265, 0.239]
    cylinder_discrete_score_without_obstacle =[0.842, 0.674, 0.635, 0.559, 0.533, 0.497, 0.478, 0.462, 0.442, 0.433, 0.430, 0.417]
    
    cylinder_implicit_score_with_obstacle=[0.759, 0.764, 0.808, 0.871, 0.899, 0.931, 0.953, 0.966, 0.971, 0.975, 0.975, 0.978]
    
    cylinder_implicit_score_without_obstacle=[0.766, 0.773, 0.825, 0.870, 0.902, 0.939, 0.960, 0.970, 0.973, 0.975, 0.976, 0.980]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_discrete_score_with_obstacle = [0.918, 0.91, 0.873, 0.806, 0.765, 0.69, 0.61, 0.545, 0.493, 0.448, 0.406, 0.359]
    vortex_street_discrete_score_without_obstacle =  [0.925, 0.926, 0.915, 0.892, 0.883, 0.871, 0.845, 0.82, 0.797, 0.767, 0.731, 0.693]
    
    
    vortex_street_implicit_score_with_obstacle=[0.930, 0.944, 0.957, 0.965, 0.971, 0.977, 0.982, 0.986, 0.986, 0.987, 0.988, 0.989]
        
        
    vortex_street_implicit_score_without_obstacle = [0.943, 0.957, 0.963, 0.968, 0.973, 0.978, 0.982, 0.986, 0.986, 0.987, 0.988, 0.989]
    
    label1=r"w/ obstacle"
    label2=r"w/o obstacle"
    
    
    draw_function(axes[0], xlabel_contour, xlabel_title_contour, cylinder_discrete_score_without_obstacle,cylinder_discrete_score_with_obstacle, 'lower right', (1.0, 0.5),rotation_angle,0,label1,label2, 1.1)
    
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
    
    
    cylinder_discrete_score_with_obstacle = [0.904, 0.933, 0.936, 0.962, 0.966, 0.979, 0.980, 0.980, 0.981, 0.983, 0.982, 0.985]
    cylinder_discrete_score_without_obstacle =[0.898, 0.928, 0.932, 0.959, 0.963, 0.977, 0.978, 0.979, 0.980, 0.982, 0.981, 0.985]
    
    cylinder_implicit_score_with_obstacle= [0.809, 0.894, 0.909, 0.940, 0.953, 0.972, 0.978, 0.984, 0.985, 0.986, 0.987, 0.988]
    
    cylinder_implicit_score_without_obstacle= [0.803, 0.890, 0.905, 0.940, 0.952, 0.970, 0.978, 0.984, 0.984, 0.985, 0.986, 0.987]
    
    xlabel_title_contour = r"Step size $s$"
    # ylabel_title_contour = r"\#CC"
    # Epsilon data

    # plt.rcParams.update({'font.size': 25})
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(17, 5),gridspec_kw={'width_ratios': [8, 8]})  # 3 subplots in a row
 # Relative widths of the subplots

    
    vortex_street_discrete_score_with_obstacle = [0.945, 0.911, 0.957, 0.947, 0.971, 0.976, 0.974, 0.977, 0.978, 0.981, 0.981, 0.981]
    vortex_street_discrete_score_without_obstacle =  [0.944, 0.910, 0.957, 0.947, 0.971, 0.976, 0.974, 0.977, 0.978, 0.981, 0.980, 0.980]
    
    
    vortex_street_implicit_score_with_obstacle=[0.697, 0.846, 0.893, 0.938, 0.951, 0.968, 0.975, 0.981, 0.983, 0.985, 0.987, 0.988]
        
        
    vortex_street_implicit_score_without_obstacle =[0.695, 0.845, 0.892, 0.937, 0.950, 0.968, 0.975, 0.981, 0.983, 0.985, 0.987, 0.988]
    
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
    
    
def gradient_norm_fig():
    
    rotation_angle=40
     # Contour data
    step_sizes = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960]
    xlabel_contour = [rf"$r/{s}$" for s in step_sizes]
    
    
    vortex_inr_mean = [1.70e-11, 1.62e-11, 1.60e-11, 1.55e-11, 1.58e-11, 1.64e-11, 1.64e-11, 1.79e-11, 1.85e-11, 2.04e-11, 2.21e-11, 2.15e-11]

    vortex_inr_max = [9.98e-11, 1.00e-10, 9.99e-11, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10]

    cylinder_inr_mean = [2.09e-11, 2.17e-11, 2.16e-11, 2.11e-11, 2.17e-11, 2.22e-11, 2.25e-11, 2.21e-11, 2.24e-11, 2.27e-11, 2.28e-11, 2.30e-11]

    cylinder_inr_max = [1.00e-10, 9.99e-11, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10]

    vortex_mfa_mean = [6.20e-12, 5.84e-12, 6.30e-12, 6.23e-12, 4.59e-12, 1.67e-12, 2.12e-12, 6.11e-12, 1.46e-11, 2.54e-11, 2.81e-11, 2.71e-11]

    vortex_mfa_max = [9.96e-11, 1.00e-10, 9.98e-11, 9.99e-11, 9.99e-11, 9.98e-11, 9.98e-11, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10]

    cylinder_mfa_mean = [4.19e-12, 5.42e-12, 9.83e-12, 1.61e-11, 1.74e-11, 1.38e-11, 1.28e-11, 1.38e-11, 1.64e-11, 1.94e-11, 2.06e-11, 2.06e-11]

    cylinder_mfa_max = [9.79e-11, 9.96e-11, 9.99e-11, 1.00e-10, 1.00e-10, 9.99e-11, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10]
    
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
    
    
def closed_form_fig():
    
    rotation_angle=60
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
    fig, axes = plt.subplots(1, 3, figsize=(17, 4),gridspec_kw={'width_ratios': [5, 5, 5]})  # 3 subplots in a row
    draw_function2(axes[0], xlabel_contour, xlabel_title_contour, potential_cc,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[1], xlabel_contour, xlabel_title_contour, rotation_cc,'lower right', (1, 0.1), rotation_angle,0)
    draw_function2(axes[2], xlabel_contour, xlabel_title_contour, potential_cc,'lower right', (1, 0.1), rotation_angle,0)
    fig.tight_layout(pad=0.1)
    fig.savefig(folder + 'closed_form_step_size.png', bbox_inches='tight', pad_inches=0)
    fig.savefig(folder + 'closd_form_step_size.pdf', bbox_inches='tight', pad_inches=0)
    fig.clear()


# mfa_fig()
# inr_fig()
closed_form_fig()
# inr_matching_score_fig() 

# mfa_matching_score_fig()

# gradient_norm_fig()