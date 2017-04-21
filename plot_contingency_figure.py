import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy
import parse_file
from numpy.random import shuffle
from math import ceil, exp,fabs
from scipy.special import digamma
from scipy.interpolate import interp1d
import sys
from stats_utils import *
import mutation_spectrum_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
import CBcm
import figure_utils
import stats_utils

if len(sys.argv) > 1:
    level = sys.argv[1]
else:
    level = 'gene'

populations = parse_file.complete_nonmutator_lines
FDR=0.05

##############################################################################
#
# Set up figures
#
##############################################################################

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

####################
#
# Fig 6: signatures of historical contingency
#
####################

fig6 = plt.figure(figsize=(7.2, 4.5))

outer_grid  = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)

lower_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[0.9,1],hspace=0.3, subplot_spec=outer_grid[1])

matrix_grid = gridspec.GridSpecFromSubplotSpec(1, 4,
                width_ratios=[1,1,1,0.05],
                wspace=0.1,
                subplot_spec=lower_grid[1])

#####
#
# Gene time distribution panel
#
#####

time_axis = plt.Subplot(fig6, outer_grid[0])
fig6.add_subplot(time_axis)
time_axis.set_ylabel('Appearance time')

time_axis.get_xaxis().tick_top()
time_axis.get_yaxis().tick_left()
time_axis.get_yaxis().set_tick_params(direction='out')
time_axis.get_xaxis().set_tick_params(direction='out')

time_axis.set_ylim([-2000,62000])
time_axis.set_yticks(figure_utils.time_xticks)
time_axis.set_yticklabels(figure_utils.time_xticklabels)

time_vmin = 0
time_vmax = 1
time_cmap = CBcm.CB2cm['redblue']
time_cNorm  = colors.Normalize(vmin=time_vmin, vmax=time_vmax)
time_scalarMap = cmx.ScalarMappable(norm=time_cNorm, cmap=time_cmap)

#####
#
# Mutation spectrum panel
#
#####

probability_axis = plt.Subplot(fig6, lower_grid[0])
fig6.add_subplot(probability_axis)
probability_axis.set_ylabel('Relative weight (%)')

probability_axis.spines['top'].set_visible(False)
probability_axis.spines['right'].set_visible(False)
probability_axis.get_xaxis().tick_bottom()
probability_axis.get_yaxis().tick_left()
probability_axis.get_yaxis().set_tick_params(direction='out')


#probability_axis.set_yticks(numpy.arange(0,5))

######
#
# Dispersion matrix panels
#
######
all_matrix_axis = plt.Subplot(fig6, matrix_grid[0])
fig6.add_subplot(all_matrix_axis)

all_matrix_axis.set_yticks([0.5,1.5,2.5,3.5,4.5])
all_matrix_axis.set_yticklabels(['6+','5','4','3','2'])
all_matrix_axis.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5])
all_matrix_axis.set_xticklabels(['1','2','3','4','5','6'])
all_matrix_axis.set_xlabel('# populations with mutation')
all_matrix_axis.set_ylabel('# mutations in gene')

all_matrix_axis.spines['top'].set_visible(False)
all_matrix_axis.spines['right'].set_visible(False)
all_matrix_axis.get_xaxis().tick_bottom()
all_matrix_axis.get_yaxis().tick_left()
all_matrix_axis.tick_params(axis='both', which='both',length=0)

all_matrix_axis.set_xlim([0,6.05])
all_matrix_axis.set_ylim([0,5.05])


early_matrix_axis = plt.Subplot(fig6, matrix_grid[1])
fig6.add_subplot(early_matrix_axis)

early_matrix_axis.set_yticks([0.5,1.5,2.5,3.5,4.5])
early_matrix_axis.set_yticklabels([])
early_matrix_axis.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5])
early_matrix_axis.set_xticklabels(['1','2','3','4','5','6'])
early_matrix_axis.set_xlabel('# populations with mutation')

early_matrix_axis.spines['top'].set_visible(False)
early_matrix_axis.spines['right'].set_visible(False)
early_matrix_axis.get_xaxis().tick_bottom()
early_matrix_axis.get_yaxis().tick_left()
early_matrix_axis.tick_params(axis='both', which='both',length=0)

early_matrix_axis.set_xlim([0,6.05])
early_matrix_axis.set_ylim([0,5.05])

late_matrix_axis = plt.Subplot(fig6, matrix_grid[2])
fig6.add_subplot(late_matrix_axis)

late_matrix_axis.set_yticks([0.5,1.5,2.5,3.5,4.5])
late_matrix_axis.set_yticklabels([])
late_matrix_axis.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5])
late_matrix_axis.set_xticklabels(['1','2','3','4','5','6'])
late_matrix_axis.set_xlabel('# populations with mutation')


late_matrix_axis.spines['top'].set_visible(False)
late_matrix_axis.spines['right'].set_visible(False)
late_matrix_axis.get_xaxis().tick_bottom()
late_matrix_axis.get_yaxis().tick_left()
late_matrix_axis.tick_params(axis='both', which='both',length=0)

late_matrix_axis.set_xlim([0,6.05])
late_matrix_axis.set_ylim([0,5.05])


cax = plt.Subplot(fig6, matrix_grid[3])
fig6.add_subplot(cax)

matrix_vmin = -20
matrix_vmax = 20
#cmap = CBcm.CB2cm['redblue']
matrix_cmap = pylab.get_cmap('RdBu_r')
matrix_cNorm  = colors.Normalize(vmin=matrix_vmin, vmax=matrix_vmax)
matrix_scalarMap = cmx.ScalarMappable(norm=matrix_cNorm, cmap=matrix_cmap)

####################################
#
# Supplemental Fig: Pooled distribution of appearance times
#
####################################

pooled_fig = plt.figure(figsize=(3, 1.7))

pooled_grid = gridspec.GridSpec(1, 1)

pooled_time_axis = plt.Subplot(pooled_fig, pooled_grid[0])
pooled_fig.add_subplot(pooled_time_axis)

pooled_time_axis.spines['top'].set_visible(False)
pooled_time_axis.spines['right'].set_visible(False)
pooled_time_axis.get_xaxis().tick_bottom()
pooled_time_axis.get_yaxis().tick_left()

pooled_time_axis.set_ylabel('Fraction mutations $\geq t$',fontsize=6)
pooled_time_axis.set_xlabel('Appearance time, $t$',fontsize=6)
pooled_time_axis.set_xticks(figure_utils.time_xticks)
pooled_time_axis.set_xticklabels(figure_utils.time_xticklabels)
pooled_time_axis.set_xlim([0,60000])



####################################
#
# Supplemental Fig: Temporal LRT as function of time
#
####################################

LRT_fig = plt.figure(figsize=(3, 1.7))

LRT_grid = gridspec.GridSpec(1, 1)

LRT_axis = plt.Subplot(LRT_fig, LRT_grid[0])
LRT_fig.add_subplot(LRT_axis)

LRT_axis.spines['top'].set_visible(False)
LRT_axis.spines['right'].set_visible(False)
LRT_axis.get_xaxis().tick_bottom()
LRT_axis.get_yaxis().tick_left()

LRT_axis.set_ylabel('LRT statistic, $\\Delta \ell$',fontsize=6)
LRT_axis.set_xlabel('Partition time, $t^*$',fontsize=6)
LRT_axis.set_xticks(figure_utils.time_xticks)
LRT_axis.set_xticklabels(figure_utils.time_xticklabels)
LRT_axis.set_xlim([0,55000])

####################################
#
# Supplemental Fig: 2-hit stuff
#
####################################

twohit_fig = plt.figure(figsize=(3, 1.7))

twohit_grid = gridspec.GridSpec(1, 1)

twohit_axis = plt.Subplot(twohit_fig, twohit_grid[0])
twohit_fig.add_subplot(twohit_axis)

twohit_axis.spines['top'].set_visible(False)
twohit_axis.spines['right'].set_visible(False)
twohit_axis.get_xaxis().tick_bottom()
twohit_axis.get_yaxis().tick_left()

twohit_axis.set_ylabel('Fraction of genes $\geq \Delta t$',fontsize=6)
twohit_axis.set_xlabel('Time difference, $\Delta t$', fontsize=6)
twohit_axis.set_xticks(figure_utils.time_xticks)
twohit_axis.set_xticklabels(figure_utils.time_xticklabels)
twohit_axis.set_xlim([0,60000])
twohit_axis.set_ylim([0,1])

####################################
#
# Supplemental Fig: Net missed opportunities as function of time
#
####################################

missed_opportunity_fig = plt.figure(figsize=(3, 1.7))

missed_opportunity_grid = gridspec.GridSpec(1, 1)

early_color = '#a50f15' 
late_color = '#045a8d'

missed_opportunity_axis = plt.Subplot(missed_opportunity_fig, missed_opportunity_grid[0])
missed_opportunity_fig.add_subplot(missed_opportunity_axis)

missed_opportunity_axis.spines['top'].set_visible(False)
missed_opportunity_axis.spines['right'].set_visible(False)
missed_opportunity_axis.get_xaxis().tick_bottom()
missed_opportunity_axis.get_yaxis().tick_left()

missed_opportunity_axis.set_ylabel('Net missed opportunities',fontsize=6)
missed_opportunity_axis.set_xlabel('Partition time, $t^*$',fontsize=6)
missed_opportunity_axis.set_xticks(figure_utils.time_xticks)
missed_opportunity_axis.set_xticklabels(figure_utils.time_xticklabels)
missed_opportunity_axis.set_xlim([0,55000])
#missed_opportunity_axis.set_ylim([-20,35])


#######################################
#
# Now do calculations and plot figures
#
#######################################

tstars = numpy.arange(0,111)*500

# Load convergence matrix
convergence_matrix = parse_file.parse_convergence_matrix(parse_file.data_directory+('%s_convergence_matrix.txt' % level))

# Load significant genes
parallel_genes = parse_file.parse_parallel_genes(parse_file.data_directory+('parallel_%ss.txt' % level))

# Calculate gene parallelism statistics
gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations)

# Calculate gene name, pop, and time vectors 

# All genes
all_gene_names = []
all_pops = []
all_times = []
for gene_name in convergence_matrix.keys():  
    for population in populations:
        for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
            all_gene_names.append(gene_name)
            all_pops.append(population)
            all_times.append(t)

# Multi-hit genes (ni>=2)
multihit_gene_names = []
multihit_pops = []
multihit_times = []

for gene_name in convergence_matrix.keys():  
    if gene_parallelism_statistics[gene_name]['observed']>1.5:
        for population in populations:
            for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                multihit_gene_names.append(gene_name)
                multihit_pops.append(population)
                multihit_times.append(t)

# One-hit genes (ni=1)
onehit_gene_names = []
onehit_pops = []
onehit_times = []

for gene_name in convergence_matrix.keys():  
    if fabs(gene_parallelism_statistics[gene_name]['observed']-1) < 0.5:
        for population in populations:
            for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                onehit_gene_names.append(gene_name)
                onehit_pops.append(population)
                onehit_times.append(t)

# Two-hit genes (ni=2)
twohit_gene_names = []
twohit_pops = []
twohit_times = []

for gene_name in convergence_matrix.keys():  
    if fabs(gene_parallelism_statistics[gene_name]['observed']-2) < 0.5:
        for population in populations:
            for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                twohit_gene_names.append(gene_name)
                twohit_pops.append(population)
                twohit_times.append(t)

# Significantly parallel genes
parallel_gene_names = []
parallel_pops = []
parallel_times = []

# Calculate gene name, pop, and time vectors 
for gene_name in parallel_genes:  
    for population in populations:
        for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
            parallel_gene_names.append(gene_name)
            parallel_pops.append(population)
            parallel_times.append(t)

parallel_gene_names = numpy.array(parallel_gene_names)
parallel_pops = numpy.array(parallel_pops)
parallel_times = numpy.array(parallel_times)

# Calculate total number of mutations and total number of significant mutations
ntot = len(all_gene_names)
nsig = len(parallel_gene_names)

# Helper function for converting vectors into time distribution
def calculate_time_distribution(gene_names, times):
    
    time_distribution = {}
    for gene_name, t in zip(gene_names,times):
        if gene_name not in time_distribution:
            time_distribution[gene_name] = []
        
        time_distribution[gene_name].append(t)
        
    for gene_name in time_distribution.keys():
        time_distribution[gene_name] = numpy.array(time_distribution[gene_name])
        time_distribution[gene_name].sort()
        
    return time_distribution

# Helper function for converting vectors into population distribution
def calculate_population_distribution(gene_names, pops):
    
    pop_distribution = {}
    for gene_name, pop in zip(gene_names,pops):
        if gene_name not in pop_distribution:
            pop_distribution[gene_name] = {population: 0 for population in populations}

        pop_distribution[gene_name][pop] += 1
        
    for gene_name in pop_distribution.keys():
        pop_distribution[gene_name] = numpy.array([pop_distribution[gene_name][pop] for pop in populations])
        
    return pop_distribution
    
# Helper function for converting vectors into distributions 
def calculate_early_late_time_distributions(gene_names, times, tstar=-1):
    
    time_distribution = {}
    for gene_name, t in zip(gene_names,times):
        if gene_name not in time_distribution:
            time_distribution[gene_name] = {'early': [], 'late' : []}
        
        if t<=tstar:
            time_distribution[gene_name]['early'].append(t)
        else:
            time_distribution[gene_name]['late'].append(t)
            
    return time_distribution


# Calculate pooled time distribution for all genes
pooled_all_time_distribution = numpy.array(all_times, copy=True)
pooled_all_time_distribution.sort()

# Calculate pooled time distribution for multihit genes
pooled_multihit_time_distribution = numpy.array(multihit_times, copy=True)
pooled_multihit_time_distribution.sort()

# Calculate pooled time distribution for 1-hit genes
pooled_onehit_time_distribution = numpy.array(onehit_times, copy=True)
pooled_onehit_time_distribution.sort()

# Calculate pooled time distribution for 2-hit genes
pooled_twohit_time_distribution = numpy.array(twohit_times, copy=True)
pooled_twohit_time_distribution.sort()
# Calculate gene-specific time distribution for 2-hit genes
twohit_time_distribution = calculate_time_distribution(twohit_gene_names, twohit_times)

# Calculate pooled time distribution for significantly parallel genes
pooled_parallel_time_distribution = numpy.array(parallel_times, copy=True)
pooled_parallel_time_distribution.sort()
# Calculate gene-specific time distribution for significantly parallel genes
parallel_time_distribution = calculate_time_distribution(parallel_gene_names, parallel_times)

# Calculate KS distance between gene-specific distributions and pooled subset
parallel_kss = numpy.array([stats_utils.calculate_ks_distance(parallel_time_distribution[gene_name], pooled_parallel_time_distribution) for gene_name in parallel_genes])

# Re-sort parallel genes by descending (num_hits, ks_statistic) 
parallel_genes, parallel_kss = (numpy.array(x) for x in zip(*sorted(zip(parallel_genes, parallel_kss), key=lambda pair: (gene_parallelism_statistics[pair[0]]['observed'],pair[1]), reverse=True)))

# Calculate median times of gene-specific time distributions
parallel_median_times = numpy.array([numpy.median(parallel_time_distribution[gene_name]) for gene_name in parallel_genes])

######################
#
# Plot pooled time distributions
#
######################

# use same color map as for distribution of appearance times by var-type
colors = [figure_utils.get_var_type_color(var_type) for var_type in parse_file.var_types]

all_ts, all_survivals = calculate_unnormalized_survival_from_vector(pooled_all_time_distribution, min_x=0, max_x=60000, min_p=1e-10)

pooled_time_axis.step(all_ts, all_survivals*1.0/all_survivals[0],'-',color='k',label='All genes')

onehit_ts, onehit_survivals = calculate_unnormalized_survival_from_vector(pooled_onehit_time_distribution, min_x=0, max_x=60000, min_p=1e-10)

pooled_time_axis.step(onehit_ts, onehit_survivals*1.0/onehit_survivals[0],'-',color=colors[0],label='1-hit genes', alpha=0.5)

twohit_ts, twohit_survivals = calculate_unnormalized_survival_from_vector(pooled_twohit_time_distribution, min_x=0, max_x=60000, min_p=1e-10)

pooled_time_axis.step(twohit_ts, twohit_survivals*1.0/twohit_survivals[0],'-',color=colors[1],label='2-hit genes', alpha=0.5)

multihit_ts, multihit_survivals = calculate_unnormalized_survival_from_vector(pooled_multihit_time_distribution, min_x=0, max_x=60000, min_p=1e-10)

pooled_time_axis.step(multihit_ts, multihit_survivals*1.0/multihit_survivals[0],'-',color=colors[2],label='Multi-hit genes', alpha=0.5)

parallel_ts, parallel_survivals = calculate_unnormalized_survival_from_vector(pooled_parallel_time_distribution, min_x=0, max_x=60000, min_p=1e-10)

interpolated_time_CDF = interp1d(parallel_ts, 1.0-parallel_survivals/parallel_survivals[0],kind='linear',bounds_error=True)

pooled_time_axis.step(parallel_ts, parallel_survivals*1.0/parallel_survivals[0],'-',color=colors[3],label='Significant genes',alpha=0.5)

pooled_time_axis.legend(loc='upper right', frameon=False)

###############
#
# Two-hit gene calculation
#
###############
       
# Calculate time difference between earliest and latest genes in 2-hit genes
observed_twohit_differences = numpy.array([twohit_time_distribution[gene_name][-1]-twohit_time_distribution[gene_name][0] for gene_name in twohit_gene_names])  
observed_twohit_differences.sort()

sys.stderr.write("Bootstrap resampling 2-hit genes...\t")
num_bootstraps=10000
bootstrapped_twohit_differences = []
for bootstrap_idx in xrange(0,num_bootstraps):

    bootstrapped_twohit_gene_names = numpy.array([gene_name for gene_name in twohit_gene_names])
    shuffle(bootstrapped_twohit_gene_names)
    
    bootstrapped_twohit_time_distribution = calculate_time_distribution(bootstrapped_twohit_gene_names, twohit_times)
    bootstrapped_twohit_differences.append( numpy.array([bootstrapped_twohit_time_distribution[gene_name][-1]-bootstrapped_twohit_time_distribution[gene_name][0] for gene_name in bootstrapped_twohit_gene_names]) )    

bootstrapped_twohit_differences = numpy.array(bootstrapped_twohit_differences)

null_twohit_differences = bootstrapped_twohit_differences.flatten()

# Calculate p-value for mean of distribution
observed_mean_twohit_difference = observed_twohit_differences.mean()
bootstrapped_mean_twohit_differences = bootstrapped_twohit_differences.mean(axis=1)
twohit_difference_pvalue = stats_utils.calculate_empirical_pvalue(-observed_mean_twohit_difference, -bootstrapped_mean_twohit_differences)

sys.stderr.write("Done!\n")
sys.stdout.write("Observed two-hit dt = %g, expected = %g, pvalue = %g (%d bootstraps)\n" % (observed_mean_twohit_difference, bootstrapped_mean_twohit_differences.mean(), twohit_difference_pvalue, num_bootstraps))

###############
#
# Plot two-hit gene calculation
#
###############

null_dts, null_survivals = calculate_unnormalized_survival_from_vector(null_twohit_differences, min_x=0, max_x=60000, min_p=1e-10)

observed_dts, observed_survivals = calculate_unnormalized_survival_from_vector(observed_twohit_differences, min_x=0, max_x=60000, min_p=1e-10)

twohit_axis.step(observed_dts, observed_survivals*1.0/observed_survivals[0],'-',color=parse_file.nonmutator_group_color,label='Observed')

twohit_axis.step(null_dts, null_survivals*1.0/null_survivals[0],'-',color='0.7',linewidth=0.5,label='Expected')

twohit_axis.legend(loc='upper right',frameon=False)

################
#
# Parallel gene calculation
#
################

observed_LRTs = []
for tstar in tstars:
    time_distribution = calculate_early_late_time_distributions(parallel_gene_names, parallel_times, tstar)
    early_ns = numpy.array([len(time_distribution[gene_name]['early'])*1.0 for gene_name in parallel_genes])
    late_ns = numpy.array([len(time_distribution[gene_name]['late'])*1.0 for gene_name in parallel_genes])
    observed_LRTs.append( mutation_spectrum_utils.calculate_LRT_statistic(early_ns, late_ns) )

observed_LRTs = numpy.array(observed_LRTs)
    
sys.stderr.write("Bootstrap resampling parallel genes...\n")
num_bootstraps=10000

bootstrapped_LRTs = []
bootstrapped_kss = []

for bootstrap_idx in xrange(1,num_bootstraps+1):

    if bootstrap_idx%1000==0:
        sys.stderr.write("%dk\n" % (bootstrap_idx/1000))

    bootstrapped_gene_names = numpy.array([gene_name for gene_name in parallel_gene_names])
    shuffle(bootstrapped_gene_names)
    
    bootstrapped_time_distribution = calculate_time_distribution(bootstrapped_gene_names, parallel_times)
    
    bootstrapped_kss.append( [stats_utils.calculate_ks_distance(bootstrapped_time_distribution[gene_name], pooled_parallel_time_distribution) for gene_name in parallel_genes] )     
    
    LRTs = []
    for tstar in tstars:

        bootstrapped_time_distribution = calculate_early_late_time_distributions(bootstrapped_gene_names, parallel_times, tstar)
    
        early_ns = numpy.array([len(bootstrapped_time_distribution[gene_name]['early'])*1.0 for gene_name in parallel_genes])
    
        late_ns = numpy.array([len(bootstrapped_time_distribution[gene_name]['late'])*1.0 for gene_name in parallel_genes])
    
        LRTs.append( mutation_spectrum_utils.calculate_LRT_statistic(early_ns, late_ns) )
    
    bootstrapped_LRTs.append(LRTs)

bootstrapped_kss = numpy.array(bootstrapped_kss) # bootstrap x genes matrix    
bootstrapped_LRTs = numpy.array(bootstrapped_LRTs) # bootstrap x tstars matrix

sys.stderr.write("Done!\n")

# Calculate ks pvalues
parallel_ks_pvalues = numpy.array([stats_utils.calculate_empirical_pvalue(parallel_kss[i],bootstrapped_kss[:,i]) for i in xrange(0,len(parallel_genes))])
# Calculate ks qvalues
parallel_ks_qvalues = stats_utils.calculate_qvalues(parallel_ks_pvalues)
# Calculate set of individually significant genes

sys.stdout.write("Temporal nonuniformity of significantly parallel genes (%d bootstraps):\n" % num_bootstraps)
for gene_name, qvalue in zip(parallel_genes, parallel_ks_qvalues):
    sys.stdout.write("%s: q=%g\n" % (gene_name, qvalue))

individually_significant_genes = set(parallel_genes[parallel_ks_qvalues<FDR])
# And remaining genes
nonsignificant_genes = parallel_genes[parallel_ks_qvalues>=FDR]

# Recalculate time distribution for nonsignificant genes
nonsignificant_time_distribution = {gene_name: parallel_time_distribution[gene_name] for gene_name in nonsignificant_genes}

nonsignificant_gene_names = []
nonsignificant_times = []
for gene_name in nonsignificant_genes:        
    nonsignificant_gene_names.extend( [gene_name]*len(nonsignificant_time_distribution[gene_name]) )
    nonsignificant_times.extend( nonsignificant_time_distribution[gene_name] )
    
pooled_nonsignificant_time_distribution = numpy.array(nonsignificant_times, copy=True)
pooled_nonsignificant_time_distribution.sort()

nonsignificant_kss = numpy.array([stats_utils.calculate_ks_distance( nonsignificant_time_distribution[gene_name], pooled_nonsignificant_time_distribution) for gene_name in nonsignificant_genes])

sys.stderr.write("Bootstrap resampling non-individually significant genes...\n")
bootstrapped_nonsignificant_kss = []
for bootstrap_idx in xrange(1,num_bootstraps+1):

    if bootstrap_idx%1000==0:
        sys.stderr.write("%dk\n" % (bootstrap_idx/1000))

    bootstrapped_gene_names = numpy.array([gene_name for gene_name in nonsignificant_gene_names])
    shuffle(bootstrapped_gene_names)
    
    bootstrapped_time_distribution = calculate_time_distribution(bootstrapped_gene_names, nonsignificant_times)
    
    bootstrapped_nonsignificant_kss.append( [stats_utils.calculate_ks_distance(bootstrapped_time_distribution[gene_name], pooled_nonsignificant_time_distribution) for gene_name in nonsignificant_genes] )     

bootstrapped_nonsignificant_kss = numpy.array(bootstrapped_nonsignificant_kss) # bootstrap x genes matrix    

sys.stderr.write("Done!\n")

# Calculate pvalue for global ks sum
observed_total_ks = nonsignificant_kss.sum()
bootstrapped_total_kss = bootstrapped_nonsignificant_kss.sum(axis=1)

total_ks_pvalue = stats_utils.calculate_empirical_pvalue(observed_total_ks, bootstrapped_total_kss)

sys.stdout.write("Individually significant_genes: %s\n" % (", ".join(individually_significant_genes)))

sys.stdout.write("Remaining total KS distance = %g, expected = %g (+/- %g), pvalue = %g\n" % (observed_total_ks, bootstrapped_total_kss.mean(), bootstrapped_total_kss.std(), total_ks_pvalue))

######################
#
# Plot early-late LRT as function of time
#
######################

upper_null_LRTs = []
for i in xrange(0,len(tstars)):

    LRTs = numpy.array(bootstrapped_LRTs[:,i],copy=True)
    LRTs.sort()
    
    upper_null_LRTs.append( LRTs[long(len(LRTs)*0.95)] )       

upper_null_LRTs = numpy.array(upper_null_LRTs)

LRT_axis.plot(tstars,observed_LRTs,'-',color=parse_file.nonmutator_group_color)
LRT_axis.fill_between(tstars, numpy.zeros_like(tstars), upper_null_LRTs,color='0.7')
LRT_axis.plot(tstars,upper_null_LRTs,'-',linewidth=0.25, color='0.6')
LRT_axis.set_ylim([0,observed_LRTs.max()*1.1])

observed_max_LRT = observed_LRTs.max()
bootstrapped_max_LRTs = bootstrapped_LRTs.max(axis=1)

sys.stdout.write("Max LRT at %d: Observed = %g, Expected = %g +/- %g, p=%g\n" % (tstars[observed_LRTs.argmax()], observed_max_LRT, bootstrapped_max_LRTs.mean(), bootstrapped_max_LRTs.std(), stats_utils.calculate_empirical_pvalue(observed_max_LRT, bootstrapped_max_LRTs)))

######################
#
# Plot gene-specific time distribution
#
######################

positions = []
current_position = 0
previous_block_position = 0
current_num_hits = gene_parallelism_statistics[parallel_genes[0]]['observed']
grey = True
for i in xrange(0,len(parallel_genes)):

    if gene_parallelism_statistics[parallel_genes[i]]['observed'] != current_num_hits:
        # reached the end of a block
        current_position+=2
        
        if grey:
            time_axis.fill_between( [previous_block_position-0.5,current_position-0.5],[-2000,-2000],[62000,62000],facecolor='0.85',linewidth=0.0)
            
        grey= not grey
        previous_block_position = current_position
        current_num_hits = gene_parallelism_statistics[parallel_genes[i]]['observed']

    current_position += 1
    
    positions.append(current_position)
    

gene_labels = []

for i in xrange(0,len(parallel_genes)):
    
    colorVal = time_scalarMap.to_rgba(interpolated_time_CDF(parallel_median_times[i]))
    
    time_axis.plot([positions[i]]*len(parallel_time_distribution[parallel_genes[i]]), parallel_time_distribution[parallel_genes[i]], 'o-',color=colorVal,markersize=1.5,linewidth=0.25,markeredgewidth=0)
    
    time_axis.plot([positions[i]],[parallel_median_times[i]],'k_',alpha=0.5,markersize=2)
    
    if parallel_genes[i] in individually_significant_genes:
        significance_string = "*"
    else:
        significance_string = ""
    
    gene_labels.append('%s%s (%d)' % (significance_string, parallel_genes[i], len(parallel_time_distribution[parallel_genes[i]])))

time_axis.set_xticks(positions)
time_axis.set_xticklabels(gene_labels, rotation='vertical',fontsize=3)
time_axis.set_xlim([positions[0]-1.5,positions[-1]+1])
probability_axis.set_xticks([])
probability_axis.set_xlim([positions[0]-1.5,positions[-1]+1])

#######
#
# Plot mutation spectrum for significant genes
#
#######

tstar = numpy.median(parallel_times)

early_late_time_distribution = calculate_early_late_time_distributions(parallel_gene_names, parallel_times, tstar)

all_ns = numpy.array([len(parallel_time_distribution[gene_name]) for gene_name in parallel_genes])
    
early_ns = numpy.array([len(early_late_time_distribution[gene_name]['early'])*1.0 for gene_name in parallel_genes])
    
late_ns = numpy.array([len(early_late_time_distribution[gene_name]['late'])*1.0 for gene_name in parallel_genes])

all_ps = all_ns*1.0/all_ns.sum()*nsig*1.0/ntot
early_ps = early_ns*1.0/early_ns.sum()*nsig*1.0/ntot
late_ps = late_ns*1.0/late_ns.sum()*nsig*1.0/ntot

probability_axis.plot(positions, all_ps*100, 'k.-',label='All',linewidth=0.5,markersize=3.0)
probability_axis.plot(positions, early_ps*100, '.-', color=time_scalarMap.to_rgba(interpolated_time_CDF(10000)), label='$\leq \mathrm{Median}(t)$',alpha=0.5,linewidth=0.5,markersize=3.0)
probability_axis.plot(positions, late_ps*100, '.-', color=time_scalarMap.to_rgba(interpolated_time_CDF(50000)),label='$>\mathrm{Median}(t)$',alpha=0.5,linewidth=0.5,markersize=3.0)

sys.stdout.write("Median time of significant genes: %g\n" % tstar)

probability_axis.legend(loc='upper right', frameon=False, ncol=3)

#######
#
# Missed opportunity calculation
#
#######

observed_time_distribution = calculate_time_distribution(all_gene_names, all_times)

observed_population_distribution = calculate_population_distribution(all_gene_names, all_pops)

desired_genes = sorted(observed_time_distribution.keys())

observed_population_matrix = numpy.array([observed_population_distribution[gene_name] for gene_name in desired_genes])

observed_median_times = numpy.array([numpy.median(observed_time_distribution[gene_name]) for gene_name in desired_genes])

observed_missed_opportunities = mutation_spectrum_utils.calculate_scaled_missed_opportunities_from_matrix(observed_population_matrix)

sys.stdout.write("Missed opportunities for >= 4-hit genes:")
for i in xrange(0,len(desired_genes)):
    if observed_population_matrix[i,:].sum() > 3.5:
        sys.stdout.write("%s n=%g m=%g\n" % (desired_genes[i], observed_population_matrix[i,:].sum(), observed_missed_opportunities[i]))


sys.stderr.write("Calculating excess missed opportunities as function of tstar...\t")
tstars = numpy.arange(0,111)*500
early_dms = []
late_dms = []
bootstrapped_early_dms = []
bootstrapped_late_dms = []
num_bootstraps = 10000
for tstar in tstars:
    
    # First look at genes with median time <= tstar
    early_population_matrix = observed_population_matrix[observed_median_times<=tstar,:]
    
    if early_population_matrix.shape[0]>0:
    
        early_total_m = mutation_spectrum_utils.calculate_scaled_missed_opportunities_from_matrix(early_population_matrix).sum()
        bootstrapped_early_total_ms = []
        for bootstrap_idx in xrange(1,num_bootstraps+1):
            bootstrapped_early_population_matrix = mutation_spectrum_utils.resample_population_matrix(early_population_matrix)
            bootstrapped_early_total_ms.append( mutation_spectrum_utils.calculate_scaled_missed_opportunities_from_matrix(bootstrapped_early_population_matrix).sum() )
        bootstrapped_early_total_ms = numpy.array(bootstrapped_early_total_ms)
        
        expected_early_total_m = bootstrapped_early_total_ms.mean()
        
        early_dms.append( early_total_m - expected_early_total_m )
        bootstrapped_early_dms.append( bootstrapped_early_total_ms - expected_early_total_m )

    else:
    
        early_dms.append( 0 )     
        bootstrapped_early_dms.append( numpy.zeros(num_bootstraps)*1.0 )
    
    # Now look at genes with median time > tstar
    late_population_matrix = observed_population_matrix[observed_median_times>tstar,:]
    
    if late_population_matrix.shape[0]>0:
    
        late_total_m = mutation_spectrum_utils.calculate_scaled_missed_opportunities_from_matrix(late_population_matrix).sum()
        bootstrapped_late_total_ms = []
        for bootstrap_idx in xrange(1,num_bootstraps+1):
            bootstrapped_late_population_matrix = mutation_spectrum_utils.resample_population_matrix(late_population_matrix)
            bootstrapped_late_total_ms.append( mutation_spectrum_utils.calculate_scaled_missed_opportunities_from_matrix(bootstrapped_late_population_matrix).sum() )
        bootstrapped_late_total_ms = numpy.array(bootstrapped_late_total_ms)
        
        expected_late_total_m = bootstrapped_late_total_ms.mean()
        
        late_dms.append( late_total_m - expected_late_total_m )
        bootstrapped_late_dms.append( bootstrapped_late_total_ms - expected_late_total_m ) 

    else:
    
        late_dms.append(0)     
        bootstrapped_late_dms.append( numpy.zeros(num_bootstraps)*1.0 )
    
    
early_dms = numpy.array(early_dms)
late_dms = numpy.array(late_dms)

bootstrapped_early_dms = numpy.transpose(numpy.array(bootstrapped_early_dms))
bootstrapped_late_dms = numpy.transpose(numpy.array(bootstrapped_late_dms))

sys.stderr.write("Done!\n")

# Calculate lower 95% CI on the (negative) early excess multiplicity
lower_early_dms = []
for i in xrange(0,bootstrapped_early_dms.shape[1]):
    
    null_dms = numpy.array(bootstrapped_early_dms[:,i],copy=True)
    null_dms.sort()
    lower_early_dms.append( null_dms[len(null_dms)*0.05] )

lower_early_dms = numpy.array(lower_early_dms)    
    
# Calculate upper 95% CI on the (positive) late excess multiplicity
upper_late_dms = []
for i in xrange(0, bootstrapped_late_dms.shape[1]):
    
    null_dms = numpy.array(bootstrapped_late_dms[:,i],copy=True)
    null_dms.sort()
    upper_late_dms.append( null_dms[long(len(null_dms)*0.95)] )

upper_late_dms = numpy.array(upper_late_dms)    

lower_dm = numpy.fmin(early_dms, lower_early_dms).min()*1.1
upper_dm = numpy.fmax(late_dms, upper_late_dms).max()*1.1

ddms = late_dms-early_dms
tstar_idx = ddms.argmax()
tstar = tstars[tstar_idx]
max_ddm = ddms[tstar_idx]
total_ddm = ddms[0]

bootstrapped_ddms = bootstrapped_late_dms-bootstrapped_early_dms
bootstrapped_max_ddms = bootstrapped_ddms.max(axis=1)
 
max_ddm_pvalue = stats_utils.calculate_empirical_pvalue(max_ddm, bootstrapped_max_ddms)
total_ddm_pvalue = stats_utils.calculate_empirical_pvalue(total_ddm, bootstrapped_ddms[:,0])

missed_opportunity_axis.plot(tstars,numpy.zeros_like(tstars),'k-',linewidth=0.25)
missed_opportunity_axis.fill_between(tstars, lower_early_dms, numpy.zeros_like(lower_early_dms),color=early_color,alpha=0.25)
#missed_opportunity_axis.plot(tstars, lower_early_dms, '-',color=early_color,linewidth=0.25)
missed_opportunity_axis.plot(tstars, early_dms, '-', color=early_color,label='$\leq t^*$')
missed_opportunity_axis.fill_between(tstars, numpy.zeros_like(upper_late_dms), upper_late_dms,color=late_color,alpha=0.25)
#missed_opportunity_axis.plot(tstars, upper_late_dms, '-',color=late_color,linewidth=0.25)
missed_opportunity_axis.plot(tstars, late_dms, '-', color=late_color,label='$> t^*$')

#missed_opportunity_axis.plot([tstar,tstar],[lower_dm,upper_dm],'k-',linewidth=0.25)
#missed_opportunity_axis.set_ylim([lower_dm,upper_dm])
missed_opportunity_axis.set_ylim([-15,20])
   
missed_opportunity_axis.legend(loc='upper right', frameon=False)
   
early_ymax = early_dms[tstar_idx]
ymax = late_dms[tstar_idx]
y0 = late_dms[0]

#print tstar, late_dms[tstar_idx], early_dms[tstar_idx], max_ddm, bootstrapped_max_ddms.mean(), bootstrapped_max_ddms.std(), max_ddm_pvalue 

sys.stdout.write("Total ddM = %g (p=%g, %d bootstraps)\n" % (total_ddm, total_ddm_pvalue, num_bootstraps))
sys.stdout.write("t^* for largest excess opportunities = %g (p=%g, %d bootstraps)\n" % (tstar, max_ddm_pvalue, num_bootstraps))
sys.stdout.write("Early=%g, Late=%g\n" % (early_dms[tstar_idx], late_dms[tstar_idx]))

########
#
# Now plot dispersion matrices for (0,60000), (0,tstar), and (tstar, 60000)
#
########
for tmin, tmax, matrix_axis in zip([-1000, -1000, tstar],[65000, tstar, 65000], [all_matrix_axis, early_matrix_axis, late_matrix_axis]):

    sys.stderr.write("Calculating dispersion matrix for %d<=t<=%d\t" % (tmin,tmax))

    population_matrix = observed_population_matrix[(observed_median_times>tmin)*(observed_median_times<=tmax),:]
    
    observed_num_pops = (population_matrix>0.5).sum(axis=1)
    observed_num_muts = (population_matrix).sum(axis=1)
    
    bootstrapped_num_pops = []
    # num hits is conserved in bootstraping
    for bootstrap_idx in xrange(1,num_bootstraps+1):
        
        bootstrapped_population_matrix = mutation_spectrum_utils.resample_population_matrix(population_matrix)
        
        bootstrapped_num_pops.append( (bootstrapped_population_matrix>0.5).sum(axis=1) )
    
    bootstrapped_num_pops = numpy.array(bootstrapped_num_pops)
    
    raw_count_matrix = []    # the number of counts in that tile
    excess_probability_matrix = [] # the deviation in expected probability of that tile
    missed_opportunity_matrix = [] # the number of missed opportunities in that tile

    for hlower, h, hupper in [(1.5,2,2.5),(2.5,3,3.5),(3.5,4,4.5),(4.5,5,5.5),(5.5,6,100.5)]:

        # The genes that belong to this row...
        desired_idxs = numpy.nonzero((observed_num_muts>hlower)*(observed_num_muts<hupper))[0]
        
        observed_counts = numpy.zeros(6)*1.0
        null_counts = numpy.zeros(6)*1.0
        
        for idx in desired_idxs:
        
            observed_counts[observed_num_pops[idx]-1]+=1
            bootstrapped_counts = numpy.array([(bootstrapped_num_pops[:,idx]==k).sum() for k in xrange(1,6+1)])
            null_counts += bootstrapped_counts*1.0/(bootstrapped_counts.sum())
            
        missed_opportunities = h-numpy.arange(1,h+1)
            
        raw_count_matrix.append(observed_counts)
        excess_probability_matrix.append((observed_counts-null_counts)/(null_counts.sum())*100)
        missed_opportunity_matrix.append(missed_opportunities)
        
    
    raw_count_matrix = numpy.array([row for row in reversed(raw_count_matrix)])
    excess_probability_matrix = numpy.array([row for row in reversed(excess_probability_matrix)])
    missed_opportunity_matrix = numpy.array([row for row in reversed(missed_opportunity_matrix)])

    # plot matrix
    m = matrix_axis.pcolor(excess_probability_matrix, cmap=matrix_cmap, vmin=matrix_vmin, vmax=matrix_vmax)

    # plot border
    matrix_axis.plot([0,2,2,3,3,4,4,5,5,6,6], [5,5,4,4,3,3,2,2,1,1,0], 'k-')
    # plot grid
    matrix_axis.plot([1,1],[0,5],'k-',linewidth=0.25)
    matrix_axis.plot([2,2],[0,4],'k-',linewidth=0.25)
    matrix_axis.plot([3,3],[0,3],'k-',linewidth=0.25)
    matrix_axis.plot([4,4],[0,2],'k-',linewidth=0.25)
    matrix_axis.plot([5,5],[0,1],'k-',linewidth=0.25)
    matrix_axis.plot([0,2],[4,4],'k-',linewidth=0.25)
    matrix_axis.plot([0,3],[3,3],'k-',linewidth=0.25)
    matrix_axis.plot([0,4],[2,2],'k-',linewidth=0.25)
    matrix_axis.plot([0,5],[1,1],'k-',linewidth=0.25)

    # plot missed opportunity # guides
    for row_idx in xrange(0,len(missed_opportunity_matrix)):
        for col_idx in xrange(0,len(missed_opportunity_matrix[row_idx])):
            
            missed_opportunities = missed_opportunity_matrix[row_idx][col_idx]
            raw_count = raw_count_matrix[row_idx][col_idx]
            
            x = col_idx+0.25
            y = row_idx+0.35
            
            #if matrix_axis==all_matrix_axis:
                # plot missed opportunity weight
            #    matrix_axis.text(x,y,'$m=%d$' % missed_opportunities,size=4)
            #else:
            #    matrix_axis.text(x-0.05,y,'$n_g = %d$' % raw_count,size=4)

            matrix_axis.text(x-0.05,y,'$n_g=%d$' % raw_count,size=4)


    sys.stderr.write("Done!\n")

# Plot some text
all_matrix_axis.text(0.05,5.2,'All genes ($\\Delta m = %0.1f$)' % y0)
early_matrix_axis.text(0.05,5.2,'Median time $<t^*$ ($\\Delta m_< = %0.1f$)' % early_ymax)
late_matrix_axis.text(0.05,5.2,'Median time $\geq t^*$ ($\\Delta m_> = %0.1f$)' % ymax)
#all_matrix_axis.text(5.9,4,'$m =$ missed opportunities\n $n_g =$ # genes',fontsize=4,horizontalalignment='right')
all_matrix_axis.text(5.9,4,'$n_g =$ # genes',fontsize=4,horizontalalignment='right')


matrix_scalarMap.set_array(excess_probability_matrix)
cbar = plt.colorbar(matrix_scalarMap, cax=cax,orientation='vertical',ticks=[-20,-10,0,10,20])
cbar.set_label('Excess probability, $\Delta P$ (%)',rotation=270,labelpad=10) 

####
#
# Figure labels
#
####
 
time_axis.text(-3, 70000, figure_utils.get_panel_label('a'),fontsize=6,fontweight='bold')
probability_axis.text(1, 2.3, figure_utils.get_panel_label('b'), fontsize=6, fontweight='bold')
all_matrix_axis.text(5.5,5.2,figure_utils.get_panel_label('c'), fontsize=6, fontweight='bold')
early_matrix_axis.text(5.5,5.2,figure_utils.get_panel_label('d'), fontsize=6, fontweight='bold')
late_matrix_axis.text(5.5,5.2,figure_utils.get_panel_label('e'), fontsize=6, fontweight='bold')
 
#######
#
# Save figures to disk
#
#######
if level=='gene':

    fig6.savefig( parse_file.figure_directory+"fig6.pdf",bbox_inches='tight')
    pooled_fig.savefig( parse_file.figure_directory+"supplemental_pooled_times.pdf", bbox_inches='tight')
    LRT_fig.savefig( parse_file.figure_directory+"supplemental_temporal_LRT.pdf",bbox_inches='tight')
    twohit_fig.savefig( parse_file.figure_directory+"supplemental_twohit_times.pdf",bbox_inches='tight')
    missed_opportunity_fig.savefig( parse_file.figure_directory+"supplemental_temporal_missed_opportunities.pdf",bbox_inches='tight')   
    
elif level=='operon':

    fig6.savefig( parse_file.figure_directory+"supplemental_operon_fig6.pdf",bbox_inches='tight')         
    missed_opportunity_fig.savefig( parse_file.figure_directory+"supplemental_operon_temporal_missed_opportunities.pdf", bbox_inches='tight')
    
else:
    pass