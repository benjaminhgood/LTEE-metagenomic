import numpy
import sys
import parse_file
import mutation_spectrum_utils
from math import log10,log,exp
import stats_utils
import pylab

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FDR = 0.05
nmin = 3
populations=parse_file.complete_nonmutator_lines

if len(sys.argv) > 1:
    level = sys.argv[1]
else:
    level = 'gene'

########
#
# Set up figure S10
#
########

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

fig = plt.figure(figsize=(2.7, 1.7))

grid = gridspec.GridSpec(1, 1)

pvalue_axis = plt.Subplot(fig, grid[0])
fig.add_subplot(pvalue_axis)

pvalue_axis.spines['top'].set_visible(False)
pvalue_axis.spines['right'].set_visible(False)
pvalue_axis.get_xaxis().tick_bottom()
pvalue_axis.get_yaxis().tick_left()

pvalue_axis.set_ylabel('Number of %ss' % level)
pvalue_axis.set_xlabel('$- \\log_{10} P$')
pvalue_axis.set_ylim([1,1000])
pvalue_axis.set_xlim([-1,20])

sys.stderr.write("Analyzing %s level parallelism...\n" % level)

# Load convergence matrix 
convergence_matrix = parse_file.parse_convergence_matrix(parse_file.data_directory+("%s_convergence_matrix.txt" % level))

# Calculate basic parallellism statistics
gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations)

# Calculate G score for entire gene (G=n*g)
gene_G_scores = mutation_spectrum_utils.calculate_G_scores(gene_parallelism_statistics)
pooled_G_scores = numpy.array(gene_G_scores.values(),copy=True)
pooled_G_scores.sort()
    
null_G_survival = mutation_spectrum_utils.NullGeneGSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics)
    
observed_Gs, observed_G_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_G_scores)
    
# Do same thing for multiplicity statistic
pooled_multiplicities = numpy.array([gene_parallelism_statistics[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics.keys()])
pooled_multiplicities.sort()
    
null_multiplicity_survival = mutation_spectrum_utils.NullGeneMultiplicitySurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )
    
observed_ms, observed_multiplicity_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_multiplicities)
 
# Do same thing for num hits
pooled_hits = numpy.array([gene_parallelism_statistics[gene_name]['observed'] for gene_name in gene_parallelism_statistics.keys()])
pooled_hits.sort()    
    
null_uniform_hit_survival = mutation_spectrum_utils.NullUniformGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )
    
null_hit_survival = mutation_spectrum_utils.NullGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )
    
observed_ns, observed_hit_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_hits)
 
# Do same thing for pvalues
gene_qvalues, gene_pvalues = mutation_spectrum_utils.calculate_parallelism_qvalues(gene_parallelism_statistics)
  
gene_logpvalues = mutation_spectrum_utils.calculate_parallelism_logpvalues(gene_parallelism_statistics)

pooled_pvalues = []
for gene_name in gene_logpvalues.keys():
    if gene_parallelism_statistics[gene_name]['observed']>=nmin:
        pooled_pvalues.append( gene_logpvalues[gene_name] )
pooled_pvalues = numpy.array(pooled_pvalues)
pooled_pvalues.sort()
    
    
null_pvalue_survival = mutation_spectrum_utils.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics, nmin=nmin)
    
observed_ps, observed_pvalue_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_pvalues, min_x=-4)
 
# Print counts of gene hits, just for reference
for n in xrange(0,21):
    sys.stdout.write('%d %d-hit %ss\n' % ((numpy.fabs(pooled_hits-n)<0.5).sum(), n, level))

 
 
# Tolerate 5% of things typically having an error
# G statistic version
#threshold_idx = numpy.nonzero((null_G_survival(observed_Gs)*1.0/observed_G_survival)<FDR)[0][0]
#Gstar = observed_Gs[threshold_idx] # lowest value where this is true
#num_significant = observed_G_survival[threshold_idx]
# Multiplicity version
#threshold_idx = numpy.nonzero((null_multiplicity_survival(observed_ms)*1.0/observed_multiplicity_survival)<FDR)[0][0]
#mstar = observed_ms[threshold_idx] # lowest value where this is true
#num_significant = observed_multiplicity_survival[threshold_idx]

# Num hits version
#threshold_idx = numpy.nonzero((null_hit_survival(observed_ns)*1.0/observed_hit_survival)<FDR)[0][0]
#nstar = observed_ns[threshold_idx] # lowest value where this is true
#num_significant = observed_hit_survival[threshold_idx]

# Pvalue version
threshold_idx = numpy.nonzero((null_pvalue_survival(observed_ps)*1.0/observed_pvalue_survival)<FDR)[0][0]
pstar = observed_ps[threshold_idx] # lowest value where this is true
num_significant = observed_pvalue_survival[threshold_idx]
    
    
sys.stdout.write("Found %d significant %ss (p* = %g)\n" % (num_significant, level, exp(-pstar)))
    
pvalue_axis.step(observed_ps/log(10), null_pvalue_survival(observed_ps),'-',label='Expected',color='k')

pvalue_axis.step(observed_ps/log(10), observed_pvalue_survival,'b-',label='Observed')
pvalue_axis.semilogy([pstar/log(10), pstar/log(10)],[5e-02,num_significant],'k-',linewidth=0.5)
pvalue_axis.semilogy([-3,pstar/log(10)],[num_significant, num_significant],'k-',linewidth=0.5)
pvalue_axis.semilogy([pstar/log(10)],[num_significant],'r.')

pvalue_axis.legend(loc='upper right',frameon=False)

fig.savefig(parse_file.figure_directory+"supplemental_%s_parallelism_pvalue.pdf" % level,bbox_inches='tight')
    
    
ntot = 0
nsignificant = 0
Ltot = 0
Lsignificant = 0

nonsignificant_genes = []
    
output_filename = parse_file.data_directory+("parallel_%ss.txt" % level)
output_file = open(output_filename,"w")

# print header
output_file.write(", ".join(["Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))

sys.stdout.write("-log p^* = %g\n" % pstar)
sys.stderr.write("Nonsignificant genes:\n")
    
for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):
    
    ntot += gene_parallelism_statistics[gene_name]['observed']
    Ltot += gene_parallelism_statistics[gene_name]['length']
    
    if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
    #if gene_G_scores[gene_name]>=Gstar:
    #if gene_qvalues[gene_name]<FDR:
        
        nsignificant += gene_parallelism_statistics[gene_name]['observed']
        Lsignificant += gene_parallelism_statistics[gene_name]['length']
    
        output_file.write("\n")
        output_file.write("%s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))
         
    else:
        
        if gene_parallelism_statistics[gene_name]['observed']>2:
                sys.stderr.write("%s, %0.1f, %d, %0.2f, %0.2f, %g\n" % (gene_name, gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))
             
        
        nonsignificant_genes.append(gene_name)

output_file.close()

observed_G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics, nonsignificant_genes)

sys.stdout.write("Significant %ss:\n" % level)
sys.stdout.write("n = %g (%0.2f of total), L = %g (%0.2f of total)\n" % (nsignificant, nsignificant*1.0/ntot, Lsignificant, Lsignificant*1.0/Ltot))
sys.stdout.write("Remaining total parallelism = %g (p=%g)\n" % (observed_G, pvalue))
