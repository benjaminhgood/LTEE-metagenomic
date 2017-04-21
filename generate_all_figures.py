#############################
#
# This script generates all figures for the manuscript
#
#############################

import os
import sys
import parse_file
from datetime import date

os.system('mkdir -p figures')

sys.stdout.write("Figure generation output for LTEE metagenomic sequencing project\n")
sys.stdout.write("Date: %s\n\n" % str(date.today()))

sys.stdout.write("Plotting Fig 1...\n")
sys.stdout.write("[Main trafic plot]\n")
return_val = os.system('python plot_full_trafic_plot_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Fig 1!\n\n")
else:
    sys.stdout.write("Error w/ Fig 1!\n\n")

sys.stdout.write("Plotting Fig 2 & supplements...\n") 
sys.stdout.write("[Rate of mutation accumulation vs time]\n")
return_val = os.system('python plot_rate_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Figs 2, etc!\n\n")
else:
    sys.stdout.write("Error w/ Figs 2, etc!\n\n")

sys.stdout.write("Plotting Fig 3...\n")
sys.stdout.write("[Quasi-stable clades]\n")
return_val = os.system('python plot_clade_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Fig 3!\n\n")
else:
    sys.stdout.write("Error w/ Fig 3!\n\n")

sys.stdout.write("Plotting Ara+1 3-way...\n")
return_val = os.system('python plot_timecourse.py p1_threeway.settings.txt')
if return_val==0:
    sys.stdout.write("Done w/ Ara+1 3-way!\n\n")
else:
    sys.stdout.write("Error w/ Ara+1 3-way!\n\n")

sys.stdout.write("Plotting Fig 4...\n") 
sys.stdout.write("[Within-clade dynamics]\n")
return_val = os.system('python plot_within_clade_dynamics_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Fig 4!\n\n")
else:
    sys.stdout.write("Error w/ Fig 4!\n\n")

sys.stdout.write("Plotting Supplemental Major/Minor comparison...\n") 
return_val = os.system('python plot_majority_minority_comparison.py')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental Major/Minor comparison\n\n")
else:
    sys.stdout.write("Error w/ Supplemental Major/Minor comparison\n\n")

sys.stdout.write("Plotting Supplemental dNdS...\n")
return_val = os.system('python plot_dNdS.py')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental dNdS!\n")
else:
    sys.stdout.write("Error w/ w/ Supplemental dNdS!\n")

sys.stdout.write("Plotting Supplemental vartype times\n")
return_val = os.system('python plot_vartype_times.py')
if return_val==0:
    sys.stdout.write("Done w/ var_type times!\n\n")
else:
    sys.stdout.write("Error w/ var_type times!\n\n")

sys.stdout.write("Plotting Supplemental site multiplicity\n")
return_val = os.system('python plot_allele_multiplicity.py')
if return_val==0:
    sys.stdout.write("Done w/ site multiplicity!\n\n")
else:
    sys.stdout.write("Error w/ site multiplicity!\n\n")

sys.stdout.write("Calculating parallelism matrices...\n")
return_val = os.system('python calculate_convergence_matrices.py')
if return_val==0:
    sys.stdout.write("Done w/ parallelism matrices!\n\n")
else:
    sys.stdout.write("Error w/ parallelism matrices!\n\n")

sys.stdout.write("Plotting Fig 5 & supplements...\t")
sys.stdout.write("[Genetic parallelism]\n")
return_val = os.system('python plot_parallelism_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Fig 5, etc!\n\n")
else:
    sys.stdout.write("Error w/ Fig 5, etc!\n\n")

sys.stdout.write("Calculating significantly parallel genes,\n")
sys.stdout.write("plotting Supplemental parallelism pvalue...\n") 
return_val = os.system('python calculate_parallel_genes_table.py')
if return_val==0:
    sys.stdout.write('Done w/ supplemental parallelism pvalue!\n\n')
else:
    sys.stdout.write("Error w/ supplemental parallelism pvalue!\n\n")

sys.stdout.write("Plotting Supplemental temporal multiplicity...\n")
return_val = os.system('python plot_temporal_multiplicity.py')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental temporal multiplicity!\n\n")
else:
    sys.stdout.write("Error w/ Supplemental temporal multiplicity!\n\n")

sys.stdout.write("Plotting Fig 6 & supplements...\n")
return_val = os.system('python plot_contingency_figure.py')
if return_val==0:
    sys.stdout.write("Done w/ Fig 6, etc.!\n\n")
else:
    sys.stdout.write("Error w/ Fig 6, etc.!\n\n")

sys.stdout.write("Plotting Supplemental iclR trajectories...\n") 
return_val = os.system('python plot_timecourse.py iclR.settings.txt')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental iclR trajectories!\n\n")
else:
    sys.stdout.write("Error w/ Supplemental iclR trajectories!\n\n")

sys.stdout.write("Plotting Supplemental hslU trajectories...\n") 
return_val = os.system('python plot_timecourse.py hslU.settings.txt')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental hslU trajectories!\n\n")
else:
    sys.stdout.write("Error w/ Supplemental hslU trajectories!\n\n")


sys.stdout.write("Plotting Supplemental atoS trajectories...\n") 
return_val = os.system('python plot_timecourse.py atoS.settings.txt')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental atoS trajectories...\n\n")
else:
    sys.stdout.write("Error w/ Supplemental atoS trajectories...\n\n")

sys.stdout.write("Plotting Supplemental argR trajectories...\n") 
sys.stdout.write("[argR trajectories]\n")
return_val = os.system('python plot_timecourse.py argR.settings.txt')
if return_val==0:
    sys.stdout.write("Done w/ Supplemental argR trajectories!\n\n")
else:
    sys.stdout.write("Error w/ Supplemental argR trajectories!\n\n")

sys.stdout.write("Plotting gene vs operon comparison...\n") 
return_val = os.system('python plot_gene_vs_operon.py operon')
if return_val==0:
    sys.stdout.write('Done w/ gene vs operon!\n\n')
else:
    sys.stdout.write("Error w/ gene vs operon!\n\n")


sys.stdout.write("Calculating significantly parallel operons,\n")
sys.stdout.write("plotting Supplemental parallelism pvalue...\n") 
return_val = os.system('python calculate_parallel_genes_table.py operon')
if return_val==0:
    sys.stdout.write('Done w/ operon supplemental parallelism pvalue!\n\n')
else:
    sys.stdout.write("Error w/ operon supplemental parallelism pvalue!\n\n")

sys.stdout.write("Plotting operon temporal multiplicity...\n")
return_val = os.system('python plot_temporal_multiplicity.py operon')
if return_val==0:
    sys.stdout.write("Done w/ operon temporal multiplicity!\n\n")
else:
    sys.stdout.write("Error w/ operon temporal multiplicity!\n\n")

sys.stdout.write("Plotting Operon Fig 6...\n")
return_val = os.system('python plot_contingency_figure.py operon')
if return_val==0:
    sys.stdout.write("Done w/ Operon Fig 6!\n\n")
else:
    sys.stdout.write("Error w/ Operon Fig 6!\n\n")

sys.stdout.write("Plotting parallel gene trajectories...\n") 
return_val = os.system('python plot_parallel_gene_timecourses.py')
if return_val==0:
    sys.stdout.write("Done w/ parallel gene timecourses!\n\n")
else:
    sys.stdout.write("Error w/ parallel gene timecourses!\n\n")

sys.stdout.write("\nDone with all figures!\n")

sys.stdout.write("\nCreating supplementary tables!\n")
return_val = os.system('python create_supplementary_tables.py')
if return_val==0:
    sys.stdout.write("Done w/ supplementary tables!\n")
else:
    sys.stdout.write("Error w/ supplementary tables!\n")

