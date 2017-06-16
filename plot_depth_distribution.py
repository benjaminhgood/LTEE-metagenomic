######
#
# This script plots the distribution of sequencing depth at different timepoints
#
#####

import sys
import numpy
import parse_file
from scipy.special import gammaln as loggamma
from math import log, exp

populations = parse_file.all_lines
   
coverages = []
nonmutator_coverages = []
    
for population in populations:
        
    sys.stderr.write("Processing %s...\t" % parse_file.get_pretty_name(population))

    # calculate mutation trajectories
    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    coverages.extend(population_avg_depths[population_avg_depths>=5])
    
    if population in parse_file.complete_nonmutator_lines:
        nonmutator_coverages.extend(population_avg_depths[population_avg_depths>=5])
    
    sys.stderr.write("Done!\n")   
        
coverages = numpy.array(coverages)
nonmutator_coverages = numpy.array(nonmutator_coverages)

print "All populations: n=%d mean=%g median=%g\n" % (len(coverages), coverages.mean(), numpy.median(coverages))
print "Nonmutator populations: n=%d mean=%g median=%g\n" % (len(nonmutator_coverages), nonmutator_coverages.mean(), numpy.median(nonmutator_coverages))
