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
from numpy.random import normal

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


fig = plt.figure(figsize=(3, 1.7))

# make three panels panels
outer_grid  = gridspec.GridSpec(1, 1)


axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(axis)

axis.set_xlabel('Mutations in operon')
axis.set_ylabel('Fraction in most-mutated gene')


axis.spines['top'].set_visible(False)
axis.spines['right'].set_visible(False)
axis.get_xaxis().tick_bottom()
axis.get_yaxis().tick_left()

axis.plot([0,20],[1,1],'k-',linewidth=0.25)

###
#
# Do calculations
#
###


populations=parse_file.complete_nonmutator_lines

tstars = numpy.arange(10,110)*500

sys.stderr.write("Loading convergence matrix...\t")
gene_convergence_matrix = parse_file.parse_convergence_matrix(parse_file.data_directory+("gene_convergence_matrix.txt"))
operon_convergence_matrix = parse_file.parse_convergence_matrix(parse_file.data_directory+("operon_convergence_matrix.txt"))

sys.stderr.write("Done!\n")

# Calculate gene parallelism statistics
gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(gene_convergence_matrix,populations)
operon_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(operon_convergence_matrix,populations)

num_plotted = 0
num_off_diagonal = 0

for operon_name in operon_parallelism_statistics.keys():
    
    genes = operon_name.split(";")
    
    if len(genes) < 1.5:
        continue
    
    n_operon = operon_parallelism_statistics[operon_name]['observed']
    
    if n_operon<1.5:
        continue
        
    num_plotted += 1
    
    
    n_genes = []
    for gene_name in operon_name.split(";"):
        if gene_name in gene_parallelism_statistics:
            n_genes.append( gene_parallelism_statistics[gene_name]['observed'] )
            
    if len(n_genes)==0:
        sys.stderr.write("Big problem!")
        sys.stdout.write("Big problem!")
     
       
    max_n_gene = max(n_genes)
    
    if max_n_gene > n_operon:
        # Shouldn't happen!
        sys.stderr.write("%s n_operon=%g max_n_gene=%g\n" % (operon_name, n_operon, max_n_gene))
    elif n_operon > 5.5 and max_n_gene < 5.5:
        num_off_diagonal+=1
        
    axis.plot([n_operon+normal(0,1)*0.1],[max_n_gene*1.0/n_operon+normal(0,1)*(1.0/n_operon)*0.1],'.',markersize=2)
axis.set_xlim([1.5,16.5])
axis.set_ylim([0,1.1])
sys.stdout.write("%d operons hit >=6 times with all genes hit <5 times\n" % num_off_diagonal)

sys.stderr.write("Saving figure...\t")
fig.savefig(parse_file.figure_directory+"supplemental_gene_vs_operon.pdf",bbox_inches='tight')
sys.stderr.write("Done!\n")
