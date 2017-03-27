#
# This script parses a combined list of primitive indel and junction evidence
# and tries to merge junctions into actual mutation events. The resulting file
# can then be fed (along with the average depth and SNP files) into the next 
# step of the pipeline (annotate_pvalues).  
#

import sys
import numpy
import pylab
import parse_file
from operator import itemgetter, attrgetter, methodcaller
from math import fabs
import bz2

input_filename = sys.argv[1]
indel_filename = sys.argv[2]

# where to read the input junctions from
file = bz2.BZ2File(input_filename,"r")

# where to write the output indels to
indel_file = bz2.BZ2File(indel_filename,"w")

# list of junctions to use to call indels
mutations = []

# how far away two JC evidences are before they are combined into one?
match_distance_threshold = 50000

for line in file:
    items = line.split(",")
    chromosome = items[0]
    location = long(items[1])
    allele = items[2].strip()
    
    if '->' in allele:
        # a snp! don't print
        continue
        
    if not allele.startswith('junction'):
        # not a junction
        # write, but don't process further
        indel_file.write(line)
        continue
    
    
    times = numpy.array([float(subitem) for subitem in items[3].split()])
    alts = numpy.array([float(subitem) for subitem in items[4].split()])
    depths = numpy.array([float(subitem) for subitem in items[5].split()])
    mutations.append((chromosome, location, allele, times, alts, depths))
file.close()

# filter mutation list
# to make sure every junction meets minimal conditions
filtered_mutations = []
for mutation in mutations:
    chromosome, location, allele, times, alts, depths = mutation
    frequencies = alts/(depths+(depths==0))
    if (alts > 1).sum() > 2 and ((alts>1)*(frequencies > 0.05)).sum() > 0:
        filtered_mutations.append(mutation)    
mutations = filtered_mutations

    
# Add fake "mutations" to to the mutation list corresponding to junctions in the reference.
# This allows us to pair junctions with the IS element on the other end in case of an IS mediated-deletion.
# We'll set the initial frequency of these mutations to 100% so we can tell that they are fake later. 
repeat_names, repeat_starts, repeat_ends, repeat_complements = parse_file.parse_repeat_list()
for i in xrange(0,len(repeat_names)):

    # make two junctions, one for start, one for end
    chromosome = mutations[0][0]
    location = repeat_starts[i]
    gene_name = repeat_names[i].split(":")[1]
    if repeat_complements[i]:
        side = 'end'
    else:
        side = 'start'
    allele = 'junction_%d_%d_%s_%s' % (repeat_starts[i],-1,gene_name,side)
    
    
    times = mutations[0][3]
    depths = numpy.ones_like(times)*1e06
    alts = numpy.ones_like(times)*1e06
    
    mutations.append((chromosome, location, allele, times, alts, depths ))
    
    
    location = repeat_ends[i]
    if repeat_complements[i]:
        side = 'start'
    else:
        side = 'end'
    allele = 'junction_%d_%d_%s_%s' % (repeat_ends[i],1,gene_name,side)
    
    mutations.append((chromosome, location, allele, times, alts, depths ))

# sort mutation list by location in genome
mutations = [mutation for mutation in sorted(mutations, key=itemgetter(1))] 

# annotate which junctions are IS mutations
junction_list = []
for mutation in mutations:
    chromosome, location, allele, times, alts, depths = mutation
    
    # parse junction string in alt_allele
    items = allele.split("_")
    if items[3].startswith("IS"): # an IS mutation
        position = long(items[1])
        strand = long(items[2])
        repeat = items[3]
        if items[4]=='start':
            side = 1
        else:
            side = -1
            
        junction_list.append((True,position,strand,repeat,side))
        
    else: # not an IS mutation
        position_1 = long(items[1])
        strand_1 = long(items[2])
        position_2 = long(items[3])
        strand_2 = long(items[4])
        junction_list.append((False,position_1,strand_1,position_2,strand_2))

# try to identify best matches between IS pairs        
match_matrix = numpy.zeros((len(junction_list),len(junction_list)))        
for i in xrange(0,len(junction_list)):
    junction_1 = junction_list[i]
    
    if junction_1[0]:
        # it involves a repeat element
        
        distances = []
        for j in xrange(0,len(junction_list)):
            junction_2 = junction_list[j]
            if junction_1[0] == junction_2[0] and junction_1[2]==-1*junction_2[2] and junction_1[3]==junction_2[3] and junction_1[4]==-1*junction_2[4]:
                distances.append(fabs(junction_2[1]-junction_1[1]))
            else:
                # if no good juction, make it infinitely far away!
                distances.append(4e09)
                
        distances = numpy.array(distances)
        
        # put a 1 in the closest entry
        # (might not be mutual)
        j = distances.argmin()
        if distances[j] < match_distance_threshold:
            #sys.stderr.write("%g\n" % distances[j])
            match_matrix[i,j] = 1
        
indel_mutations = []
matched = numpy.zeros(len(junction_list))

already_processed = [False]*len(junction_list)
                    
for i in xrange(0,len(junction_list)):

    chromosome, location, allele, times, alts, depths = mutations[i]    

    is_fake_mutation_i = ((alts[0]/(depths[0]+(depths[0]==0))) > 0.9) and (depths[0] > 1e04)
    
    
    # don't add an entry for a fake mutation
    if is_fake_mutation_i:
        continue
        
    # don't add an entry for something that has already been matched
    if already_processed[i]:
        continue

    # we're definitely going to add an entry for this one, 
    # so mark it as processed
    already_processed[i] = True

    if match_matrix[i].sum() > 0.5:
        # it has a match somewhere!
        # get it:
        j = match_matrix[i].argmax() 
        other_alts = mutations[j][4]
        other_depths = mutations[j][5]
        
        # is this one a fake mutation?
        is_fake_mutation_j = ((other_alts[0]/(other_depths[0]+(other_depths[0]==0))) > 0.9) and (other_depths[0] > 1e04)
        
        left_junction = junction_list[i]
        right_junction = junction_list[j]
        
        if is_fake_mutation_j:
            # fake mutation on the other side
            new_alts = alts
            new_depths = depths
            #indel = ('MOB', left_junction[1], right_junction[1], left_junction[3], left_junction[4])    
            #new_alt_allele = "_".join([str(item) for item in indel])
            new_alt_allele = allele
        
        elif match_matrix[j,i] > 0.5:
            # is reciprocal match! add alts and depths
            new_alts = alts+other_alts
            new_depths = depths+other_depths
            indel = ('MOB', left_junction[1], right_junction[1], left_junction[3], left_junction[4])    
            new_alt_allele = "_".join([str(item) for item in indel])
        
            # mark the new one as processed too so that we do not double count
            already_processed[j]=True 
        
        else:
            
            # not a fake one and not a reciprocal match
            # use match to call MOB event, but don't merge alts and depths
            
            new_alts = alts
            new_depths = depths
            
            #indel = ('MOB', left_junction[1], right_junction[1], left_junction[3], left_junction[4])    
            #new_alt_allele = "_".join([str(item) for item in indel])
            new_alt_allele = allele
         
            
    else:
        
        # no matches, just leave as is
        
        new_alts = alts
        new_depths = depths
        new_alt_allele = allele
        
    
    indel_mutations.append((chromosome, location, new_alt_allele, times, new_alts, new_depths))
        
            
# remove fake ones that are still hanging around 
# these should already be gone? 
final_indel_mutations = []
for mutation in indel_mutations:
    A = mutation[4][0]*1.0
    D = mutation[5][0]
    #sys.stderr.write("%g %g\n" % (A,D))
    if D < 1e05:
        final_indel_mutations.append(mutation)

for chromosome, location, allele, times, alts, depths in final_indel_mutations:
    indel_file.write(", ".join([chromosome, str(location), allele, " ".join([str(t) for t in times]), " ".join([str(a) for a in alts]), " ".join([str(d) for d in depths])]))
    indel_file.write("\n")
    
indel_file.close()
    