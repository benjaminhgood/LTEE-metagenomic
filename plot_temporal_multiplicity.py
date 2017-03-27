import sys
import numpy

from numpy.random import shuffle
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import parse_file
import mutation_spectrum_utils
import stats_utils
import figure_utils

if len(sys.argv)>1:
    level=sys.argv[1]
else:
    level='gene'

##################################
#
# Set up figures
#
##################################

mpl.rcParams['font.size'] = 6.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

all_color = 'k'
early_color = '#a50f15' # dark red
late_color = '#045a8d'

#######################################################
#
# Fig S13: Multiplicity vs time plot
#
#######################################################


fig = plt.figure(figsize=(7, 2))

# make three panels panels
outer_grid  = gridspec.GridSpec(1, 2, width_ratios=[4,3], hspace=0.2)


g_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(g_axis)

g_axis.set_xlabel('Partition time, $t^*$')
g_axis.set_ylabel('Parallelism G-score')

g_axis.set_ylim([-0.06,0.06])
g_axis.set_xticks(figure_utils.time_xticks)
g_axis.set_xticklabels(figure_utils.time_xticklabels)
g_axis.set_xlim([0,55000])

g_axis.spines['top'].set_visible(False)
g_axis.spines['right'].set_visible(False)
g_axis.get_xaxis().tick_bottom()
g_axis.get_yaxis().tick_left()


multiplicity_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(multiplicity_axis)

multiplicity_axis.set_xlabel('%s multiplicity, $m$' % level.capitalize())
multiplicity_axis.set_ylabel('Fraction mutations $\geq m$')

multiplicity_axis.set_xlim([1,1e02])
multiplicity_axis.set_ylim([1e-03,1])

multiplicity_axis.spines['top'].set_visible(False)
multiplicity_axis.spines['right'].set_visible(False)
multiplicity_axis.get_xaxis().tick_bottom()
multiplicity_axis.get_yaxis().tick_left()

###
#
# Do calculations
#
###


populations=parse_file.complete_nonmutator_lines

tstars = numpy.arange(0,110)*500

sys.stderr.write("Loading convergence matrix...\t")
convergence_matrix = parse_file.parse_convergence_matrix(parse_file.data_directory+("%s_convergence_matrix.txt" % level))
sys.stderr.write("Done!\n")


# Calculate gene parallelism statistics
gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations)

# Calculate gene scores
gene_g_scores = mutation_spectrum_utils.calculate_g_scores(gene_parallelism_statistics)

null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)
null_cdf = (lambda m: 1-null_survival_function(m))

ms = []
gs = []
ts = []

for gene_name in convergence_matrix.keys():   
    m = gene_parallelism_statistics[gene_name]['multiplicity']
    g = gene_g_scores[gene_name]
      
    for population in populations:
        for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
            ms.append(m)
            gs.append(g)
            ts.append(t)

ms = numpy.array(ms)
gs = numpy.array(gs)
ts = numpy.array(ts)

global_g = gs.mean()
global_G = gs.sum()

early_dGs = []
late_dGs = []

sys.stderr.write("Calculating fractional g scores...\t")
for tstar in tstars:
    
    
    g1s = gs[ts<=tstar]
    g2s = gs[ts>tstar]
    
    early_dGs.append((g1s-global_g).sum()/global_G)
    late_dGs.append((g2s-global_g).sum()/global_G)
    
sys.stderr.write("Done!\n")
    
early_dGs = numpy.array(early_dGs)
late_dGs = numpy.array(late_dGs)

sys.stderr.write("Bootstrap resampling...\t")

bootstrapped_ts = numpy.array(ts,copy=True)
num_bootstraps = 10000

bootstrapped_early_dGs = []
bootstrapped_late_dGs = []
for bootstrap_idx in xrange(0,num_bootstraps):

    shuffle(bootstrapped_ts)
    
    bootstrapped_early_dGs.append([(gs[bootstrapped_ts<=tstar]-global_g).sum()/global_G for tstar in tstars])
    bootstrapped_late_dGs.append([(gs[bootstrapped_ts>tstar]-global_g).sum()/global_G for tstar in tstars])
    

bootstrapped_early_dGs = numpy.array( bootstrapped_early_dGs )
bootstrapped_late_dGs = numpy.array( bootstrapped_late_dGs )
sys.stderr.write("Done!\n")

sys.stderr.write("Calculating confidence intervals...\t")
upper_early_dGs = []
lower_late_dGs = []
for i in xrange(0,len(tstars)):

    dGs = numpy.array(bootstrapped_early_dGs[:,i],copy=True)
    dGs.sort()
    
    upper_early_dGs.append( dGs[long(len(dGs)*0.975)] )
    
    dGs = numpy.array(bootstrapped_late_dGs[:,i],copy=True)
    dGs.sort()
    
    lower_late_dGs.append( dGs[long(len(dGs)*0.025)] )
    
upper_early_dGs = numpy.array(upper_early_dGs)
lower_late_dGs = numpy.array(lower_late_dGs)
sys.stderr.write("Done!\n")

bootstrapped_max_ddGs = (bootstrapped_early_dGs-bootstrapped_late_dGs).max(axis=1)
observed_max_ddG = (early_dGs-late_dGs).max()


tstar = tstars[(early_dGs-late_dGs).argmax()]

sys.stdout.write("Plotting multiplicity distribution for tstar = %g\n" % tstar)    
sys.stdout.write("ddG = %g, expected = %g +/- %g, p=%g\n" % (observed_max_ddG, bootstrapped_max_ddGs.mean(), bootstrapped_max_ddGs.std(), stats_utils.calculate_empirical_pvalue(observed_max_ddG, bootstrapped_max_ddGs)) )

###
#
# Plot figures
#
###


g_axis.plot(tstars, numpy.zeros_like(tstars),'-',color=all_color,label='All',linewidth=0.5)

g_axis.fill_between(tstars, numpy.zeros_like(tstars), upper_early_dGs, color='0.7',linewidth=0)
g_axis.plot(tstars, upper_early_dGs, color='0.6',linewidth=0.25)
g_axis.fill_between(tstars, lower_late_dGs, numpy.zeros_like(tstars), color='0.7',linewidth=0)
g_axis.plot(tstars, lower_late_dGs, color='0.6',linewidth=0.25)

g_axis.plot(tstars, early_dGs,'-',color=early_color, label='$\leq t^*$')
g_axis.plot(tstars, late_dGs, '-',color=late_color, label='$> t^*$')

g_axis.legend(loc='upper right',frameon=False)

early_survival_ms, early_survivals = stats_utils.calculate_unnormalized_survival_from_vector(ms[ts<=tstar], min_x=0.1, max_x=100)
late_survival_ms, late_survivals = stats_utils.calculate_unnormalized_survival_from_vector(ms[ts>tstar], min_x=0.1, max_x=100)
all_survival_ms, all_survivals = stats_utils.calculate_unnormalized_survival_from_vector(ms, min_x=0.1, max_x=100)

theory_ms = numpy.logspace(0,2,100)
multiplicity_axis.loglog(theory_ms, null_survival_function(theory_ms),color='0.7',linewidth=0.5)
multiplicity_axis.step(all_survival_ms, all_survivals*1.0/all_survivals[0],color=all_color,linewidth=0.5)
multiplicity_axis.step(early_survival_ms, early_survivals*1.0/early_survivals[0],color=early_color)
multiplicity_axis.step(late_survival_ms, late_survivals*1.0/late_survivals[0],color=late_color)

if level=='gene':
    fig.savefig( parse_file.figure_directory+"supplemental_temporal_multiplicity.pdf",bbox_inches='tight')
   
elif level=='operon':
    fig.savefig( parse_file.figure_directory+"supplemental_operon_temporal_multiplicity.pdf",bbox_inches='tight')
    
