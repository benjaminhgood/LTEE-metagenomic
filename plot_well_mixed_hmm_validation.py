import numpy
import pylab
import sys
from numpy.random import binomial, choice
import matplotlib as mpl
from numpy.random import binomial,normal
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import parse_file

import uniform_prior_well_mixed_hmm as well_mixed_model

import stats_utils
import timecourse_utils
import figure_utils

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


coverage = 50

def calculate_true_appearance_fixation_time(times, freqs, freq_threshold=0.05):
    
    polymorphic_idxs = (freqs>=freq_threshold)*(freqs<=(1-freq_threshold))
            
    if polymorphic_idxs.sum() == 0:
        
        appearance_time = times[freqs.argmax()]-250
        fixation_time = times[freqs.argmax()]+250
        transit_time = fixation_time-appearance_time
    
    else:
    
        # do regular calculation
        appearance_time = times[polymorphic_idxs][0]-250
        
        if polymorphic_idxs[-1]:
            # Stays polymorphic until the end
            fixation_time = 1000000 # never fixes
            transit_time = times[-1]+250 - appearance_time
        else:
            fixation_time = times[polymorphic_idxs][-1]+250
            transit_time = fixation_time - appearance_time
        
    return appearance_time, fixation_time, transit_time

def parse_timecourse(filename):
    file = open(filename,"r")
    file.readline()
    times = numpy.array([float(item) for item in file.readline().split()])
    total_snp_trajectories = []
    fixed_mutations = []
    for line in file:
        fixed_mutation_items = line.split(":")[2].split()
        for item in fixed_mutation_items:
            fixed_mutations.append(float(item.split(",")[0]))
        snp_map = {}
        line_items = line.split(":")[3].split("+")
        total_snp_data = line_items[:-1]
        for i in xrange(0, len(times)):
            t = times[i]
            snp_data = total_snp_data[i]
            for snp in snp_data.split(";")[:-1]:
                items = snp.split(",")
                label = long(items[0])
                frequency = float(items[1])  
                if not label in snp_map:
                   snp_map[label] = {}
                snp_map[label][i] = frequency
    
        snp_trajectories = []
        for snp_label in snp_map.keys():
            trajectory = []
            for i in xrange(0,len(times)):
                if i in snp_map[snp_label]:

                    # add noise                   
                    #trajectory.append( binomial(coverage,snp_map[snp_label][i])*1.0/coverage)
                    trajectory.append(snp_map[snp_label][i])
                else:
                    trajectory.append(0)
            snp_trajectories.append(numpy.array(trajectory))
        total_snp_trajectories.append(snp_trajectories)
    return times, total_snp_trajectories, numpy.array(fixed_mutations) 
    
#times, total_snp_trajectories, fixed_mutations = parse_timecourse("additional_data/macroscopic_epistasis_simulations/wiser_timecourse_beneficial_output.txt")

times, total_snp_trajectories, fixed_mutations = parse_timecourse("additional_data/macroscopic_epistasis_simulations/wiser_timecourse_output.txt")


hmm_appearance_times = []
true_appearance_times = []
hmm_transit_times = []
true_transit_times = []

focal_population_idx = 1
focal_population_empirical_freqs = None

for population_idx in xrange(0,len(total_snp_trajectories)):

    true_freqs = []
    empirical_freqs = []
    D = []
    A = []

    for freqs in total_snp_trajectories[population_idx]:
        if (freqs>0.1).any():
        
            true_freqs.append(freqs)
            
            depths = numpy.ones_like(times)*coverage*1.0
            alts = binomial(coverage, freqs)*1.0
            
            empirical_freqs.append(alts/depths)
            D.append(depths)
            A.append(alts)

    A = numpy.array(A)
    D = numpy.array(D)
        
    Pstate, Ls, p0, p1 = well_mixed_model.infer_hmm(A,D,num_iterations=10)
    
    for mutation_idx in xrange(0,len(true_freqs)):
        hmm_appearance_time, hmm_fixation_time, hmm_transit_time = timecourse_utils.calculate_appearance_fixation_time_from_hmm(times,empirical_freqs[mutation_idx],Ls[mutation_idx,:]) 
        
        true_appearance_time, true_fixation_time, true_transit_time = calculate_true_appearance_fixation_time(times, true_freqs[mutation_idx])
            
        hmm_appearance_times.append(hmm_appearance_time)
        true_appearance_times.append(true_appearance_time)
        hmm_transit_times.append(hmm_transit_time)
        true_transit_times.append(true_transit_time)
    
    if population_idx==focal_population_idx:
        focal_population_empirical_freqs = empirical_freqs
            
hmm_appearance_times = numpy.array(hmm_appearance_times)
true_appearance_times = numpy.array(true_appearance_times)
hmm_transit_times = numpy.array(hmm_transit_times)
true_transit_times = numpy.array(true_transit_times)

transit_errors = hmm_transit_times-true_transit_times
transit_error_xs, transit_error_survivals = stats_utils.calculate_unnormalized_survival_from_vector(transit_errors)


appearance_errors = hmm_appearance_times-true_appearance_times
error_xs, error_survivals = stats_utils.calculate_unnormalized_survival_from_vector(appearance_errors)

bins = numpy.arange(-20,21)*500+250


print len(total_snp_trajectories), "populations"
print len(appearance_errors), "trajectories"

####
#
# Set up figure
#
####

fig = plt.figure(figsize=(7, 3))


outer_grid = gridspec.GridSpec(2, 1, height_ratios=[0.7,1],hspace=0.4)


inner_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[1,1],
                subplot_spec=outer_grid[1], wspace=0.1) #, hspace=0.08)


trafic_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(trafic_axis)

trafic_axis.spines['top'].set_visible(False)
trafic_axis.spines['right'].set_visible(False)
trafic_axis.get_xaxis().tick_bottom()
trafic_axis.get_yaxis().tick_left()

trafic_axis.set_ylabel('Allele frequency, $f(t)$')
trafic_axis.set_xlabel('Generation, $t$')

trafic_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
trafic_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

trafic_axis.set_ylim([0,1.05])
trafic_axis.set_xlim([0,60000])

trafic_axis.set_xticks(figure_utils.time_xticks)
trafic_axis.set_xticklabels(figure_utils.time_xticklabels)

#trafic_axis.set_yticks(frequency_ticks)


appearance_axis = plt.Subplot(fig, inner_grid[0])
fig.add_subplot(appearance_axis)

appearance_axis.spines['top'].set_visible(False)
appearance_axis.spines['right'].set_visible(False)
appearance_axis.get_xaxis().tick_bottom()
appearance_axis.get_yaxis().tick_left()

appearance_axis.set_ylabel('Number of mutations')
appearance_axis.set_xlabel('Appearance time error')

appearance_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
appearance_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)


transit_axis = plt.Subplot(fig, inner_grid[1])
fig.add_subplot(transit_axis)

transit_axis.spines['top'].set_visible(False)
transit_axis.spines['right'].set_visible(False)
transit_axis.get_xaxis().tick_bottom()
transit_axis.get_yaxis().tick_left()

transit_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
transit_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)


transit_axis.set_xlabel('Transit time error')
transit_axis.set_xlim([-5000,5000])
appearance_axis.set_xlim([-5000,5000])
appearance_axis.set_ylim([0,13000])
transit_axis.set_ylim([0,13000])

transit_axis.hist(transit_errors, bins=bins, facecolor=figure_utils.nonmutator_group_color, edgecolor=figure_utils.nonmutator_group_color)
appearance_axis.hist(appearance_errors, bins=bins, facecolor=figure_utils.nonmutator_group_color, edgecolor=figure_utils.nonmutator_group_color)

for fs in focal_population_empirical_freqs:
    trafic_axis.plot(times,fs,'.-',markersize=1,linewidth=0.5)

#fig.savefig(parse_file.figure_directory+'supplemental_well_mixed_hmm_validation.pdf',bbox_inches='tight',transparent=True)
fig.savefig(parse_file.figure_directory+'supplemental_well_mixed_hmm_validation.png',bbox_inches='tight',transparent=True,dpi=300)
