import sys
import os
import parse_file

###########################################################################
#
# This script does not do any calculations. It merely copies and renames 
# existing data files to the figures folder as supplementary_table_X.csv
#
###########################################################################

sys.stderr.write("Creating Supplementary table 1...\n")
sys.stdout.write("[List of metagenomic samples]\n")
os.system('cp population_samples.csv %ssupplementary_table_1.csv' % parse_file.figure_directory)
sys.stdout.write("Done!\n")

sys.stderr.write("Creating Supplementary table 2...\n")
sys.stdout.write("[List of clonal isolates]\n")
os.system('cp clone_samples.csv %ssupplementary_table_2.csv' % parse_file.figure_directory)
sys.stdout.write("Done!\n")

sys.stderr.write("Creating Supplementary table 3...\n")
sys.stdout.write("[List of significantly parallel genes in nonmutators]\n")
os.system('cp data_files/parallel_genes.txt %ssupplementary_table_3.csv' % parse_file.figure_directory)
sys.stdout.write("Done!\n")

sys.stderr.write("Creating Supplementary table 4...\n")
sys.stdout.write("[List of significantly parallel operons in nonmutators]\n")
os.system('cp data_files/parallel_operons.txt %ssupplementary_table_4.csv' % parse_file.figure_directory)
sys.stdout.write("Done!\n")