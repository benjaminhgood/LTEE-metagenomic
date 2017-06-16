##############################################################################
#
# This script generates Fig. 2 in the main text
# as well as supplementary Figures S1, S2, S3, and S4.
# 
##############################################################################

import sys
import numpy
from math import log10
import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import parse_file
import timecourse_utils
from numpy.random import normal, choice, shuffle

import figure_utils
import stats_utils

mpl.rcParams['font.size'] = 5
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

populations = parse_file.complete_nonmutator_lines + parse_file.mutator_lines

##############################################################################
#
# Set up figures
#
##############################################################################

##############################################################################
#
# First set up Fig. 2 (rates of fitness and mutation accumulation with time)
#
##############################################################################

fig2 = plt.figure(figsize=(5, 2))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 3, width_ratios=[0.95, 0.55, 0.05], wspace=0.43)

trajectory_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.4,1],
                subplot_spec=outer_grid[0], hspace=0.05) #, hspace=0.08)

middle_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.6,1],
                subplot_spec=outer_grid[1], hspace=0.5)

###################
#
# Legend Panel
#
###################

legend_axis = plt.Subplot(fig2, outer_grid[2])
fig2.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])

#############################
#
# Fitness vs time panel
#
#############################

fitness_axis = plt.Subplot(fig2, trajectory_grid[0])
fig2.add_subplot(fitness_axis)

fitness_axis.set_ylabel('Fitness, $X(t)$ (%)',labelpad=2)
fitness_axis.set_ylim([0,45])

fitness_axis.spines['top'].set_visible(False)
fitness_axis.spines['right'].set_visible(False)
fitness_axis.get_xaxis().tick_bottom()
fitness_axis.get_yaxis().tick_left()

#############################
#
# Mutations vs time panel
#
#############################
mut_axis = plt.Subplot(fig2, trajectory_grid[1])
fig2.add_subplot(mut_axis)

mut_axis.set_xlabel('Generation, $t$',labelpad=2)
mut_axis.set_ylabel('Avg # mutations, $M(t)$',labelpad=2)

mut_axis.set_xlim([-1000, 61000])
mut_axis.set_xticks(figure_utils.time_xticks)
mut_axis.set_xticklabels(figure_utils.time_xticklabels)

fitness_axis.set_xticks(figure_utils.time_xticks)
fitness_axis.set_xticklabels(["" for label in figure_utils.time_xticklabels])
fitness_axis.set_xlim([-1000,61000])
fitness_axis.set_yticks([0,10,20,30,40])

fitness_axis.tick_params(axis='y', direction='out',length=3,pad=1)

mut_axis.tick_params(axis='y', direction='out',length=3,pad=1)
mut_axis.tick_params(axis='x', direction='out',length=3,pad=1)

mut_axis.set_ylim([0,119])

mut_axis.spines['right'].set_visible(False)
mut_axis.get_xaxis().tick_bottom()
mut_axis.get_yaxis().tick_left()

#############################
#
# Mutator mutations vs time inset
#
#############################    
  
#mutator_mut_axis = inset_axes(mut_axis, width="30%", height="30%", borderpad=0, bbox_to_anchor=(-0.60,-0.05,1, 1), bbox_transform=mut_axis.transAxes) # at top left

mutator_mut_axis = inset_axes(mut_axis, width="23%", height="23%", borderpad=0, bbox_to_anchor=(-0.01,-0.625,1, 1), bbox_transform=mut_axis.transAxes)
#-0.025

mutator_mut_axis.set_xlim([-1000,61000])
mutator_mut_axis.set_xticks([0,20000,40000,60000])
mutator_mut_axis.set_xticklabels(['0','20k','40k','60k'])
mutator_mut_axis.set_ylim([0,3000])      
mutator_mut_axis.set_yticks([0,1000,2000,3000])
mutator_mut_axis.set_yticklabels(['0','1k','2k','3k'])

mutator_mut_axis.fill_between(numpy.array([0,61000]),numpy.array([0,0]), numpy.array([3000,3000]),color='w',zorder=20)

    
cl = pylab.getp(mutator_mut_axis, 'xmajorticklabels')
pylab.setp(cl, fontsize=4) 
cl = pylab.getp(mutator_mut_axis, 'ymajorticklabels')
pylab.setp(cl, fontsize=4) 
mutator_mut_axis.get_yaxis().tick_left()
mutator_mut_axis.get_xaxis().tick_bottom()
[i.set_linewidth(0.75) for i in mutator_mut_axis.spines.itervalues()]

#############################
#
# Derivative vs time panel
#
#############################
velocity_axis = plt.Subplot(fig2, middle_grid[0])
fig2.add_subplot(velocity_axis)

velocity_axis.set_xlabel('Generation, $t$',labelpad=2,fontsize=4)
velocity_axis.set_ylabel('$dM / dt$ ($\\times 10^3$)',labelpad=2)

velocity_axis.set_xlim([-1000, 61000])
velocity_axis.set_xticks(figure_utils.time_xticks)
velocity_axis.set_xticklabels(figure_utils.time_xticklabels)
velocity_axis.set_yticks([0,1,2,3])
velocity_axis.set_ylim([0,3.5])

velocity_axis.spines['top'].set_visible(False)
velocity_axis.spines['right'].set_visible(False)
velocity_axis.get_xaxis().tick_bottom()
velocity_axis.get_yaxis().tick_left()

velocity_axis.tick_params(axis='y', direction='out',length=3,pad=1,labelsize=4)
velocity_axis.tick_params(axis='x', direction='out',length=3,pad=1,labelsize=4)

#############################
#
# Fixed vs avg mutations panel
#
#############################

fixation_axis = plt.Subplot(fig2, middle_grid[1])

fig2.add_subplot(fixation_axis)

fixation_axis.set_xlabel('Avg # mutations',labelpad=2,fontsize=4)
fixation_axis.set_ylabel('  Fixed mutations',labelpad=4,rotation=270,fontsize=4)
fixation_axis.yaxis.set_label_position("right")

fixation_axis.set_ylim([0,115])
fixation_axis.set_xlim([0,115])

fixation_axis.spines['top'].set_visible(False)
fixation_axis.spines['left'].set_visible(False)

fixation_axis.get_xaxis().tick_bottom()
fixation_axis.get_yaxis().tick_right()

fixation_axis.tick_params(axis='y', direction='out',length=3,pad=1,labelsize=4)
fixation_axis.tick_params(axis='x', direction='out',length=3,pad=1,labelsize=4)

fixation_axis.plot([0,115],[0,115],'-',linewidth=0.5)

##############################################################################
#
# Now set up Supplemental Fig (Panels A and B of Fig. 2 redone with W measure)
#
##############################################################################

w_fig = plt.figure(figsize=(2.25, 2))

# Set up grids to hold figure panels
trajectory_grid = gridspec.GridSpec(2, 1, height_ratios=[0.4,1], hspace=0.05)

#############################
#
# Fitness vs time panel
#
#############################

w_fitness_axis = plt.Subplot(w_fig, trajectory_grid[0])
w_fig.add_subplot(w_fitness_axis)

w_fitness_axis.set_ylabel('Fitness, $W(t)$',labelpad=2)
w_fitness_axis.set_ylim([1,2.15])

w_fitness_axis.spines['top'].set_visible(False)
w_fitness_axis.spines['right'].set_visible(False)
w_fitness_axis.get_xaxis().tick_bottom()
w_fitness_axis.get_yaxis().tick_left()

#############################
#
# Mutations vs time panel
#
#############################
w_mut_axis = plt.Subplot(w_fig, trajectory_grid[1])
w_fig.add_subplot(w_mut_axis)

w_mut_axis.set_xlabel('Generation, $t$',labelpad=2)
w_mut_axis.set_ylabel('Avg # mutations, $M(t)$',labelpad=2)

w_mut_axis.set_xlim([-1000, 61000])
w_mut_axis.set_xticks(figure_utils.time_xticks)
w_mut_axis.set_xticklabels(figure_utils.time_xticklabels)

w_fitness_axis.set_xticks(figure_utils.time_xticks)
w_fitness_axis.set_xticklabels(["" for label in figure_utils.time_xticklabels])
w_fitness_axis.set_xlim([-1000,61000])

w_fitness_axis.tick_params(axis='y', direction='out',length=3,pad=1)

w_mut_axis.tick_params(axis='y', direction='out',length=3,pad=1)
w_mut_axis.tick_params(axis='x', direction='out',length=3,pad=1)

w_mut_axis.set_ylim([0,119])

w_mut_axis.spines['right'].set_visible(False)
w_mut_axis.get_xaxis().tick_bottom()
w_mut_axis.get_yaxis().tick_left()

##############################################################################
#
# Now set up Supplemental Fig (avg_s vs time)
#
##############################################################################

s_fig = plt.figure(figsize=(2, 1.5))

# Set up grids to hold figure panels
trajectory_grid = gridspec.GridSpec(1, 1)

#############################
#
# Fitness / Mutations vs time panel
#
#############################
s_axis = plt.Subplot(s_fig, trajectory_grid[0])
s_fig.add_subplot(s_axis)

s_axis.set_xlabel('Generation, $t$',labelpad=2)
s_axis.set_ylabel('Fitness per mutation, $X(t)/M(t)$ (%)',labelpad=2)

s_axis.set_xlim([-1000, 61000])
s_axis.set_xticks(figure_utils.time_xticks)
s_axis.set_xticklabels(figure_utils.time_xticklabels)

s_axis.tick_params(axis='y', direction='out',length=3,pad=1)
s_axis.tick_params(axis='x', direction='out',length=3,pad=1)

s_axis.set_ylim([0,5])

s_axis.spines['top'].set_visible(False)
s_axis.spines['right'].set_visible(False)
s_axis.get_xaxis().tick_bottom()
s_axis.get_yaxis().tick_left()

##############################################################################
#
# Now set up Supplemental Figure (Fitness and mutation accumulation from 40k-60k)
#
##############################################################################

late_fig = plt.figure(figsize=(4.3, 1.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.4)

late_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.7,1], hspace=0.05,subplot_spec=outer_grid[0])

early_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.7,1], hspace=0.05,subplot_spec=outer_grid[1])

late_mut_axis = plt.Subplot(late_fig, late_grid[1])
late_fig.add_subplot(late_mut_axis)

late_mut_axis.set_xlabel('Generation, $t$',labelpad=2)
late_mut_axis.set_ylabel('Mutation gain',labelpad=2)

late_mut_axis.set_xlim([40000, 61000])
late_mut_axis.set_ylim([-1,32.5])
xticks = [10000*i for i in xrange(4,7)]
xticklabels = ['%dk' % (10*i) for i in xrange(4,7)]
late_mut_axis.set_xticks(xticks)
late_mut_axis.set_xticklabels(xticklabels)

late_mut_axis.spines['top'].set_visible(False)
late_mut_axis.spines['right'].set_visible(False)
late_mut_axis.get_xaxis().tick_bottom()
late_mut_axis.get_yaxis().tick_left()

late_fitness_axis = plt.Subplot(late_fig, late_grid[0])
late_fig.add_subplot(late_fitness_axis)

late_fitness_axis.set_ylabel('Fitness gain\n$\Delta X$ (%)',labelpad=2)

late_fitness_axis.spines['top'].set_visible(False)
late_fitness_axis.spines['right'].set_visible(False)
late_fitness_axis.get_xaxis().tick_bottom()
late_fitness_axis.get_yaxis().tick_left()
late_fitness_axis.set_xticks(xticks)
late_fitness_axis.set_xticklabels(["" for tick in xticks])
late_fitness_axis.set_yticks([0,1,2,3,4])

late_fitness_axis.set_xlim([40000,61000])
late_fitness_axis.set_ylim([0,4.5])

early_mut_axis = plt.Subplot(late_fig, early_grid[1])
late_fig.add_subplot(early_mut_axis)

early_mut_axis.set_xlabel('Generation, $t$',labelpad=2)
early_mut_axis.set_ylabel('Mutation gain',labelpad=2)

early_mut_axis.set_xlim([0, 10000])
early_mut_axis.set_ylim([-1,32.5])
xticks = [0,5000,10000]
xticklabels = ['0','5k','10k']
early_mut_axis.set_xticks(xticks)
early_mut_axis.set_xticklabels(xticklabels)

early_mut_axis.spines['top'].set_visible(False)
early_mut_axis.spines['right'].set_visible(False)
early_mut_axis.get_xaxis().tick_bottom()
early_mut_axis.get_yaxis().tick_left()

early_fitness_axis = plt.Subplot(late_fig, early_grid[0])
late_fig.add_subplot(early_fitness_axis)

early_fitness_axis.set_ylabel('Fitness gain\n$\Delta X$ (%)',labelpad=2)

early_fitness_axis.spines['top'].set_visible(False)
early_fitness_axis.spines['right'].set_visible(False)
early_fitness_axis.get_xaxis().tick_bottom()
early_fitness_axis.get_yaxis().tick_left()
early_fitness_axis.set_xticks(xticks)
early_fitness_axis.set_xticklabels(["" for tick in xticks])

#early_fitness_axis.set_yticks([0,1,2,3,4])

early_fitness_axis.set_xlim([0,10000])
early_fitness_axis.set_ylim([0,30])


##############################################################################
#
# Now set up Supplemental Fig (Temporal correlations in mutation accumulation)
#
##############################################################################

correlation_fig = plt.figure(figsize=(5, 1.5))

inner_grid_1 = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.3) 

correlation_mut_axis = plt.Subplot(correlation_fig, inner_grid_1[0])
correlation_fig.add_subplot(correlation_mut_axis)
correlation_mut_axis.set_xlabel('Generation, $t$',labelpad=2)
correlation_mut_axis.set_ylabel('Avg # of mutations, $M$',labelpad=2)

correlation_mut_axis.set_xlim([0, 61000])
correlation_mut_axis.set_xticks(figure_utils.time_xticks)
correlation_mut_axis.set_xticklabels(figure_utils.time_xticklabels)
correlation_mut_axis.set_ylim([0,115])

correlation_mut_axis.spines['top'].set_visible(False)
correlation_mut_axis.spines['right'].set_visible(False)
correlation_mut_axis.get_xaxis().tick_bottom()
correlation_mut_axis.get_yaxis().tick_left()

correlation_pvalue_axis = plt.Subplot(correlation_fig, inner_grid_1[1])
correlation_fig.add_subplot(correlation_pvalue_axis)
correlation_pvalue_axis.set_xlabel('Between line stddev, $\sigma_M$', labelpad=2)
correlation_pvalue_axis.set_ylabel('Fraction bootstraps $\geq \sigma_M$', labelpad=2)

correlation_pvalue_axis.set_ylim([1e-04,2])
correlation_pvalue_axis.set_xlim([0,25])

correlation_pvalue_axis.spines['top'].set_visible(False)
correlation_pvalue_axis.spines['right'].set_visible(False)
correlation_pvalue_axis.get_xaxis().tick_bottom()
correlation_pvalue_axis.get_yaxis().tick_left()
#correlation_pvalue_axis.set_xticks([])

##############################################################################
#
# Now set up Fig. S5 (Fixed vs avg mutations in mutators)
#
##############################################################################

mutator_fixation_fig = plt.figure(figsize=(2, 1.5))
mutator_fixation_grid = gridspec.GridSpec(1, 1)
mutator_fixation_axis = plt.Subplot(mutator_fixation_fig, mutator_fixation_grid[0])

mutator_fixation_fig.add_subplot(mutator_fixation_axis)

mutator_fixation_axis.set_xlabel('Avg # mutations ($\\times 10^3$)',labelpad=2)
mutator_fixation_axis.set_ylabel('Fixed mutations ($\\times 10^3$)',labelpad=8,rotation=270)
mutator_fixation_axis.yaxis.set_label_position("right")

mutator_fixation_axis.spines['top'].set_visible(False)
mutator_fixation_axis.spines['left'].set_visible(False)

mutator_fixation_axis.get_xaxis().tick_bottom()
mutator_fixation_axis.get_yaxis().tick_right()

mutator_fixation_axis.tick_params(axis='y', direction='out',length=3,pad=1)
mutator_fixation_axis.tick_params(axis='x', direction='out',length=3,pad=1)

mutator_fixation_axis.set_ylim([-0.025,3])
mutator_fixation_axis.set_xlim([0,3.025])      
 
mutator_fixation_axis.plot([0,3],[0,3],'k-',linewidth=0.5)

##############################################################################
#
# Now set up Supplemental Fig (distribution of transit times)
#
##############################################################################

transit_time_fig = plt.figure(figsize=(2, 1.5))
transit_time_grid = gridspec.GridSpec(1, 1)
transit_time_axis = plt.Subplot(transit_time_fig, transit_time_grid[0])

transit_time_fig.add_subplot(transit_time_axis)

transit_time_axis.set_xlabel('Transit time, $\\Delta t$')
transit_time_axis.set_ylabel('Fraction mutations $\\geq \\Delta t$')

transit_time_axis.set_xticks(figure_utils.time_xticks)
transit_time_axis.set_xticklabels(figure_utils.time_xticklabels)

transit_time_axis.spines['top'].set_visible(False)
transit_time_axis.spines['right'].set_visible(False)
transit_time_axis.get_xaxis().tick_bottom()
transit_time_axis.get_yaxis().tick_left()



##############################################################################
#
# Now perform actual calculations used for figures
#
##############################################################################

fitness_trajectories = {}
W_fitness_trajectories = {}
late_fitness_trajectories = {}

mutation_trajectories = {}
fixed_mutation_trajectories = {}

transit_times = {}

sys.stderr.write("Loading fitness data...\n")

sys.stderr.write("Loading Wiser fitness data...\t")
trajectories, line_data = parse_file.parse_ancestor_fitnesses()
for population in trajectories.keys():
    
    fitness_ts = trajectories[population][0]
    fitness_xs = trajectories[population][1]
    fitness_ws = trajectories[population][2]
    
    fitness_trajectories[population] = (numpy.array(fitness_ts), numpy.array(fitness_xs))
    W_fitness_trajectories[population] = (fitness_ts, fitness_ws)
sys.stderr.write("Done!\n")

sys.stderr.write("Loading late fitness gains...\t")
measured_populations,individual_50k_dxs, individual_50k_std_dxs, avg_50k_dx, std_50k_dx, std_50k_dx_measurement, individual_60k_dxs, individual_60k_std_dxs, avg_60k_dx, std_60k_dx, std_60k_dx_measurement = parse_file.get_all_40k_50k_fitnesses()
for i in xrange(0,len(measured_populations)):
    population = measured_populations[i]
    
    if not population in parse_file.complete_nonmutator_lines:
        continue
    
    t40 = 40000
    t50 = 50000+normal(0,500)
    t60 = 60000+normal(0,500)
    dx40 = 0
    dx50 = individual_50k_dxs[i]
    dx60 = individual_60k_dxs[i]
    std_dx40 = 0
    std_dx50 = individual_50k_std_dxs[i]
    std_dx60 = individual_60k_std_dxs[i]
    
    late_fitness_ts = numpy.array([t40, t50, t60])
    late_fitness_xs = numpy.array([dx40, dx50, dx60])
    late_fitness_std_xs = numpy.array([std_dx40, std_dx50, std_dx60])
    
    late_fitness_trajectories[population] = (late_fitness_ts, late_fitness_xs, late_fitness_std_xs)
sys.stderr.write("Done!\n")


sys.stderr.write("Loading mutation data...\n")
for population in populations:
    
    sys.stderr.write("Processing %s...\t" % parse_file.get_pretty_name(population))

    # calculate mutation trajectories
    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
        
    times = mutations[0][10]
    Ms = numpy.zeros_like(times)*1.0
    fixed_Ms = numpy.zeros_like(times)*1.0
    
    transit_times[population] = []
    
    for mutation_idx in xrange(0,len(mutations)):
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
        Ls = haplotype_trajectories[mutation_idx]
        state_Ls = state_trajectories[mutation_idx]
        
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]
        
        t0,tf,transit_time = timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
        transit_times[population].append(transit_time)
        
          
        interpolating_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
   
        fs = interpolating_function(times)
        fs[fs<0]=0
        
        # Record 
        Ms += fs
        if masked_state_Ls[-1] in parse_file.well_mixed_fixed_states:
            fixed_Ms += (times>=tf)
        

    fixed_mutation_trajectories[population] = (times, fixed_Ms)
    mutation_trajectories[population] = (times, Ms)

    sys.stderr.write("analyzed %d mutations!\n" % len(mutations))

sys.stderr.write("Done!\n")    

sys.stderr.write("Postprocessing trajectories...\t") 

# Calculate average trajectories for complete nonmutator populations
avg_Mts, avg_Ms = timecourse_utils.average_trajectories([mutation_trajectories[population] for population in parse_file.complete_nonmutator_lines])
avg_Xts, avg_Xs = timecourse_utils.average_trajectories([fitness_trajectories[population] for population in parse_file.complete_nonmutator_lines])

max_freqs = []
final_freqs = []

intercepts = {}

slope_trajectories = {}
coarse_slope_trajectories = {}

max_nonmutator_M = max([mutation_trajectories[population][1][-1] for population in parse_file.complete_nonmutator_lines])
min_nonmutator_M = min([mutation_trajectories[population][1][-1] for population in parse_file.complete_nonmutator_lines])
max_mutator_M = max([mutation_trajectories[population][1][-1] for population in parse_file.mutator_lines])
min_mutator_M = min([mutation_trajectories[population][1][-1] for population in parse_file.mutator_lines])


for population in populations:

    times, Ms = mutation_trajectories[population]
    
    # Estimate of M(t) at generation 40k
    target_idx = numpy.fabs(times-40000).argmin()
    intercepts[population]=numpy.median(Ms[target_idx-2:target_idx+3])

    dt = 2500
    slopes = []
    for t in times:
        window_idxs = (times>=t-dt)*(times<=t+dt)
        window_ts = times[window_idxs]
        window_Ms = Ms[window_idxs]
        slope =  timecourse_utils.estimate_slope(window_ts, window_Ms)
        slopes.append(slope)
    slopes = numpy.array(slopes)*1000    
    slope_trajectories[population] = (times, slopes)
    
    coarse_times = numpy.array([10000,20000,30000,40000,50000,60000])*1.0
    coarse_slopes = []
    for tmax in coarse_times:
        tmin = tmax-10000
        window_idxs = (times>=tmin)*(times<=tmax)
        window_ts = times[window_idxs]
        window_Ms = Ms[window_idxs]
        slope =  timecourse_utils.estimate_slope(window_ts, window_Ms)
        coarse_slopes.append(slope)
    coarse_slopes = numpy.array(coarse_slopes)    
    coarse_slope_trajectories[population] = (coarse_times, coarse_slopes)

# Calculate average slope in sliding windows
avg_dMts, avg_dMs = timecourse_utils.average_trajectories([slope_trajectories[population] for population in parse_file.complete_nonmutator_lines])

# Bootstrap confidence intervals for average slope
num_bootstraps = 10000
dM_map = {}
for bootstrap_idx in xrange(0, num_bootstraps):

    bootstrapped_populations = choice(parse_file.complete_nonmutator_lines, len(parse_file.complete_nonmutator_lines), replace=True)
    bootstrapped_avg_dMts, bootstrapped_avg_dMs = timecourse_utils.average_trajectories([slope_trajectories[population] for population in bootstrapped_populations])
    for t,dM in zip(bootstrapped_avg_dMts, bootstrapped_avg_dMs):
        if t not in dM_map:
            dM_map[t] = []
        dM_map[t].append(dM)
    
upper_dMs = []
lower_dMs = []
for t in avg_dMts:
    dM_map[t].sort()
    lower_dMs.append( dM_map[t][long(len(dM_map[t])*0.025)])
    upper_dMs.append( dM_map[t][long(len(dM_map[t])*0.975)])

upper_dMs = numpy.array(upper_dMs)
lower_dMs = numpy.array(lower_dMs)
sys.stderr.write("Done!\n")


# Calculate coarse-grained trajectories
sys.stderr.write("Calculate coarse-grained trajectories...\n")
coarse_mutation_trajectories = {}
dMs = []
dMs_without_p1 = []

freq_dep_idxs = []
non_freq_dep_idxs = []
freq_dep_idxs_without_p1 = []
non_freq_dep_idxs_without_p1 = []

for population in parse_file.complete_nonmutator_lines:

    coarse_dMts, coarse_dMs = coarse_slope_trajectories[population]
    coarse_dMs *= 10000
    
    coarse_Mts = numpy.hstack([[0], coarse_dMts])
    coarse_Ms = numpy.hstack([[0], coarse_dMs]).cumsum()
    coarse_mutation_trajectories[population] = (coarse_Mts, coarse_Ms)
    
    if population in set(['p2','p4']):
        non_freq_dep_idxs.append(True)
        freq_dep_idxs.append(False)
    else:
        non_freq_dep_idxs.append(False)
        freq_dep_idxs.append(True)
    
    if population != 'p1':
        dMs_without_p1.append(coarse_dMs)
        
        
        if population in set(['p2','p4']):
            non_freq_dep_idxs_without_p1.append(True)
            freq_dep_idxs_without_p1.append(False)
        else:
            non_freq_dep_idxs_without_p1.append(False)
            freq_dep_idxs_without_p1.append(True)
    
        
    dMs.append(coarse_dMs)

dMs = numpy.array(dMs)
dMs_without_p1 = numpy.array(dMs_without_p1)

freq_dep_idxs = numpy.array(freq_dep_idxs)
non_freq_dep_idxs = numpy.array(non_freq_dep_idxs)
freq_dep_idxs_without_p1 = numpy.array(freq_dep_idxs_without_p1)
non_freq_dep_idxs_without_p1 = numpy.array(non_freq_dep_idxs_without_p1)

def sample_null_dM(dMs):
    
    transposed_dMs = []
    for i in xrange(0,dMs.shape[1]):
        changes = dMs[:,i]*1.0
        shuffle(changes)
        transposed_dMs.append(changes)
    
    return numpy.transpose(numpy.array(transposed_dMs))
        
observed_stddev = (dMs.sum(axis=1)).std()
observed_stddev_without_p1 = (dMs_without_p1.sum(axis=1)).std()

observed_freq_dep_difference = (dMs.sum(axis=1))[freq_dep_idxs].mean()-(dMs.sum(axis=1))[non_freq_dep_idxs].mean()

observed_freq_dep_difference_without_p1 = (dMs_without_p1.sum(axis=1))[freq_dep_idxs_without_p1].mean()-(dMs_without_p1.sum(axis=1))[non_freq_dep_idxs_without_p1].mean()


bootstrapped_stddevs = []
bootstrapped_stddevs_without_p1 = []

bootstrapped_freq_dep_differences = []
bootstrapped_freq_dep_differences_without_p1 = []


num_bootstraps = 100000
for bootstrap_idx in xrange(2,2+num_bootstraps):
    #pylab.figure(bootstrap_idx)
    bootstrapped_dMs = sample_null_dM(dMs)
    bootstrapped_stddevs.append( (bootstrapped_dMs.sum(axis=1)).std() )
    
    bootstrapped_freq_dep_differences.append( (bootstrapped_dMs.sum(axis=1))[freq_dep_idxs].mean()-(bootstrapped_dMs.sum(axis=1))[non_freq_dep_idxs].mean() )
    
    bootstrapped_dMs_without_p1 = sample_null_dM(dMs_without_p1)
    bootstrapped_stddevs_without_p1.append( (bootstrapped_dMs_without_p1.sum(axis=1)).std() )

    bootstrapped_freq_dep_differences_without_p1.append( (bootstrapped_dMs_without_p1.sum(axis=1))[freq_dep_idxs_without_p1].mean()-(bootstrapped_dMs_without_p1.sum(axis=1))[non_freq_dep_idxs_without_p1].mean() )

bootstrapped_stddevs.sort()
bootstrapped_stddevs_without_p1.sort()
bootstrapped_freq_dep_differences.sort()
bootstrapped_freq_dep_differences_without_p1.sort()

bootstrapped_stddevs = numpy.array(bootstrapped_stddevs)
bootstrapped_stddevs_without_p1 = numpy.array(bootstrapped_stddevs_without_p1)
bootstrapped_freq_dep_differences = numpy.array(bootstrapped_freq_dep_differences)
bootstrapped_freq_dep_differences_without_p1 = numpy.array(bootstrapped_freq_dep_differences_without_p1)



pvalue = (bootstrapped_stddevs>=observed_stddev).sum()*1.0/(len(bootstrapped_stddevs))

pvalue_without_p1 = (bootstrapped_stddevs_without_p1>=observed_stddev_without_p1).sum()*1.0/(len(bootstrapped_stddevs_without_p1))

freq_dep_pvalue = stats_utils.calculate_empirical_pvalue(observed_freq_dep_difference, bootstrapped_freq_dep_differences)

freq_dep_pvalue_without_p1 = stats_utils.calculate_empirical_pvalue(observed_freq_dep_difference_without_p1, bootstrapped_freq_dep_differences_without_p1)


sys.stdout.write("Rate correlation pvalue: p=%g\n" % pvalue)
sys.stdout.write("without p1: p=%g\n" % pvalue_without_p1)

sys.stdout.write("Freq dep difference: %g (p=%g)\n" % (observed_freq_dep_difference, freq_dep_pvalue))
sys.stdout.write("without p1: %g (p=%g)\n" % (observed_freq_dep_difference_without_p1, freq_dep_pvalue_without_p1))

sys.stderr.write("Done!\n")

##############################################################################
#
# Now plot the resulting figures!
#
##############################################################################

population_idx = 0
mutator_idx = 0
nonmutator_idx = 0

for population in populations:
    
    Xts,Xs = fitness_trajectories[population]
    Wts,Ws = W_fitness_trajectories[population]
    Mts,Ms = mutation_trajectories[population]
    fixed_Mts, fixed_Ms = fixed_mutation_trajectories[population]
    dMts, dMs = slope_trajectories[population]
    intercept = intercepts[population]
    
    # calculate fitness / mutations at times where both are measured
    common_time_set = (set(Xts) & set(Mts))
    
    sts = []
    ss = []
    for t in sorted(common_time_set):
        
        x_t_idx = numpy.fabs(Xts-t).argmin()
        m_t_idx = numpy.fabs(Mts-t).argmin()
        
        sts.append(t)
        ss.append( Xs[x_t_idx]/(Ms[m_t_idx]+(Ms[m_t_idx]==0) ) )
        
        #if population in parse_file.complete_nonmutator_lines:
        #    print t, Xts[t_idx], Mts[t_idx], Xs[t_idx], Ms[t_idx], ss[-1]
    
    sts = numpy.array(sts)
    ss = numpy.array(ss)
    
    transit_times[population].sort()
    transit_times[population] = numpy.array(transit_times[population])
    
    dts, dt_survival = stats_utils.calculate_unnormalized_survival_from_vector(transit_times[population], min_x=0)

    
    if population in parse_file.complete_nonmutator_lines:
        # We're dealing with a non-mutator population    
        
        colorVal = parse_file.get_line_color(population)
        linestyle = 'o-'
        zorder = 12-nonmutator_idx
        
        nonmutator_idx += 1
        
        fixation_axis.plot(Ms, fixed_Ms, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)    
        
        late_Xts, late_Xs, late_std_Xs = late_fitness_trajectories[population]
        
        late_mut_axis.plot(Mts, Ms-intercept, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)  
        
        late_fitness_axis.plot(late_Xts, late_Xs*100, '.-', color=colorVal, markersize=3, linewidth=0.25)
        for t,X,std_X in zip(late_Xts, late_Xs, late_std_Xs)[1:]:    
            late_fitness_axis.plot([t,t], [(X-std_X)*100, (X+std_X)*100], '-', color=colorVal, linewidth=0.5)
        
        early_mut_axis.plot(Mts, Ms, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)  
        
        early_fitness_axis.plot(Xts, Xs*100, linestyle, color=colorVal, markersize=1, linewidth=0.5,zorder=zorder, markeredgewidth=0)
        
          
        coarse_Mts, coarse_Ms = coarse_mutation_trajectories[population]
        
        correlation_mut_axis.plot(Mts, Ms, '-', color='0.8',linewidth=0.25)
        correlation_mut_axis.plot(coarse_Mts, coarse_Ms, '.-', color=colorVal, markersize=3,linewidth=1)
    
        s_axis.plot(sts[1:], ss[1:]*100, linestyle, color=colorVal, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0.0)
        
    else:
        # 'mutator!'    
        colorVal = parse_file.get_line_color(population)
        zorder = 6-mutator_idx
        linestyle = 's-' 
        mutator_idx += 1 
        mutator_fixation_axis.plot(Ms/1000.0, fixed_Ms/1000.0, '-', color=colorVal, alpha=1,linewidth=0.5,zorder=zorder+50, markeredgewidth=0)         
    
    population_idx += 1
    
    fitness_axis.plot(Xts, Xs*100, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0.0) 
    w_fitness_axis.plot(Wts, Ws, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0.0) 
    mut_axis.plot(Mts, Ms, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)      
    w_mut_axis.plot(Mts, Ms, linestyle, color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)      
    
    legend_axis.plot([-2,-1],[-2,-1],linestyle,color=colorVal,alpha=1,markersize=2,linewidth=0.5,markeredgewidth=0,label=parse_file.get_pretty_name(population))
    
    # inset
    mutator_mut_axis.plot(Mts, Ms, '-', color=colorVal, alpha=1,markersize=2,linewidth=0.5,zorder=50+zorder)

    transit_time_axis.step(dts, dt_survival/dt_survival[0], color=colorVal, alpha=1, markersize=1,linewidth=0.5,zorder=zorder, markeredgewidth=0)

mut_axis.plot(avg_Mts, avg_Ms,'k-',linewidth=2,zorder=13)
mut_axis.plot(avg_Mts, avg_Ms,'w-',linewidth=1,zorder=14)
w_mut_axis.plot(avg_Mts, avg_Ms,'k-',linewidth=2,zorder=13)
w_mut_axis.plot(avg_Mts, avg_Ms,'w-',linewidth=1,zorder=14)

velocity_axis.fill_between(avg_dMts,lower_dMs,upper_dMs,facecolor='0.8',linewidth=0.3,edgecolor='0.7')
velocity_axis.plot(avg_dMts,avg_dMs,'w-',linewidth=1,zorder=14,path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()],label='$\\Delta t = 5k$')

legend_axis.legend(loc='upper center',frameon=False,fontsize=4,numpoints=1,ncol=1,handlelength=1)   

transit_time_axis.semilogy([-1,1])
transit_time_axis.set_xlim([0,60000])
transit_time_axis.set_ylim([1e-02,1])

# calculate fitness / mutations at times where both are measured
common_time_set = (set(avg_Xts) & set(avg_Mts))
    
avg_sts = []
avg_ss = []
for t in sorted(common_time_set):
        
    x_t_idx = numpy.fabs(avg_Xts-t).argmin()
    m_t_idx = numpy.fabs(avg_Mts-t).argmin()  
        
    avg_sts.append(t)
    avg_ss.append( avg_Xs[x_t_idx]/(avg_Ms[m_t_idx]+(avg_Ms[m_t_idx]==0) ) )
    
avg_sts = numpy.array(avg_sts)
avg_ss = numpy.array(avg_ss)

s_axis.plot(avg_sts[1:], avg_ss[1:]*100,'k-',linewidth=2,zorder=13)
s_axis.plot(avg_sts[1:], avg_ss[1:]*100,'w-',linewidth=1,zorder=14)

###
#
# Rate correlation pvalue
#
###

correlation_pvalue_axis.step(bootstrapped_stddevs_without_p1, 1-numpy.arange(0,len(bootstrapped_stddevs_without_p1))*1.0/len(bootstrapped_stddevs_without_p1), where='post',color='0.8')

correlation_pvalue_axis.step(bootstrapped_stddevs, 1-numpy.arange(0,len(bootstrapped_stddevs))*1.0/len(bootstrapped_stddevs), where='post',color='k')

correlation_pvalue_axis.semilogy([observed_stddev, observed_stddev],[0.1/len(bootstrapped_stddevs),pvalue],'r-o',markersize=3,linewidth=0.5,label='All nonmutator')
correlation_pvalue_axis.plot([0,observed_stddev],[pvalue,pvalue],'r-',linewidth=0.5)

correlation_pvalue_axis.semilogy([observed_stddev_without_p1, observed_stddev_without_p1],[0.1/len(bootstrapped_stddevs_without_p1), pvalue_without_p1],'r-s',markersize=3,linewidth=0.5,alpha=0.5,label='Without Ara+1')
correlation_pvalue_axis.plot([0,observed_stddev_without_p1],[pvalue_without_p1, pvalue_without_p1],'r-',linewidth=0.5,alpha=0.5)
correlation_pvalue_axis.legend(loc='upper right',frameon=False,fontsize=4,numpoints=1,ncol=1,handlelength=1)   

##
#
# Adding figure labels in Fig. 2
#
##

fitness_axis.text(1500, 38, figure_utils.get_panel_label('a'),fontsize=6,fontweight='bold')
mut_axis.text(1500, 100, figure_utils.get_panel_label('b'), fontsize=6, fontweight='bold')
velocity_axis.text(3000, 3.05, figure_utils.get_panel_label('c'), fontsize=6, fontweight='bold')
fixation_axis.text(6, 101, figure_utils.get_panel_label('d'), fontsize=6, fontweight='bold')


sys.stderr.write("Saving figures...\t")
fig2.savefig(parse_file.figure_directory+'fig2.pdf',bbox_inches='tight',transparent=True)
w_fig.savefig(parse_file.figure_directory+'supplemental_w_comparison.pdf',bbox_inches='tight',transparent=True)
s_fig.savefig(parse_file.figure_directory+'supplemental_fitness_per_mutation.pdf',bbox_inches='tight',transparent=True)
s_fig.savefig(parse_file.figure_directory+'supplemental_fitness_per_mutation.png',bbox_inches='tight',transparent=True,dpi=300)

late_fig.savefig(parse_file.figure_directory+'supplemental_late_rate_comparison.pdf',bbox_inches='tight',transparent=True)
late_fig.savefig(parse_file.figure_directory+'supplemental_late_rate_comparison.png',bbox_inches='tight',transparent=True,dpi=300)
correlation_fig.savefig(parse_file.figure_directory+'supplemental_rate_correlation.pdf',bbox_inches='tight',transparent=True)
mutator_fixation_fig.savefig(parse_file.figure_directory+'supplemental_mutator_fixation.pdf',bbox_inches='tight',transparent=True)
transit_time_fig.savefig(parse_file.figure_directory+'supplemental_well_mixed_transit_times.pdf',bbox_inches='tight',transparent=True)

sys.stderr.write("Done!\n")