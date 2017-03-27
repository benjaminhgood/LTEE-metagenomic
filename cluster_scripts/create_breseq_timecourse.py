# This script takes the output from our second breseq run 
# (the "rebreseq" step") and calculates a timecourse for 
# each junction. Because breseq will call junctions slightly
# differently at different timepoints, this script includes
# a fuzzy matching algorithm for merging junction candidates
# that are likely to be the same. 

# the idea is that allele reads get added together and references get averaged

import numpy
import sys
import population_parameters
from operator import itemgetter, attrgetter, methodcaller
from parse_file import parse_repeat_list, get_repeat_idx, get_closest_repeat_idx
from math import fabs

# File containing the reference fasta file
reference_filename = sys.argv[1]
# The population to compile timecourses for
population = sys.argv[2]
# The output (gd) files from breseq 
gd_filenames = sys.argv[3:]

# How far can two simple indels be before they are merged
INDEL_EDGE_TOLERANCE = 0
# How long can a simple indel be before it is treated differently?
INDEL_LENGTH_TOLERANCE = 100
# How far can two similar IS junctions be before they are merged
REPEAT_EDGE_TOLERANCE = 20

# How far can two similar inversion elements be before they are merged
OTHER_EDGE_TOLERANCE = 20


INDEL = 0
REPEAT = 1
OTHER = 2

# Construct the list of samples for this population
full_sample_list = population_parameters.sample_names[population]
full_sample_times = population_parameters.sample_times[population]

sample_list = []
gd_filename_map = {}
#parse sample names from filenames
for gd_filename in gd_filenames:
    sample_name = gd_filename.split("/")[2].strip()
    if sample_name.startswith(population):
        sample_name = sample_name.split("_",1)[-1]

    sample_list.append(sample_name)
    gd_filename_map[sample_name] = gd_filename

sample_times = [full_sample_times[full_sample_list.index(sample)] for sample in sample_list]

# make sure the list is sorted
sample_times, sample_list = (list(x) for x in zip(*sorted(zip(sample_times,sample_list))))

sample_times_set = sorted(set(sample_times))

# Load reference fasta 
fasta_file = open(reference_filename,"r")
fasta_file.readline() # header
reference_sequence = []
for line in fasta_file:
    reference_sequence.append(line.strip())
fasta_file.close()
reference_sequence = "".join(reference_sequence)

mutation_map = {}
repeat_data = parse_repeat_list("data/lenski_reference/REL606.6.gbk")
repeat_names, repeat_starts, repeat_ends, repeat_complements = repeat_data

mutation_str_map = {'short':{}, 'long':{}}
str_mutation_map = {'short':{}, 'long':{}}

indel_alleles = {}
  
def fuzzy_match(mutation_1, mutation_2, version='long'):
       
    if mutation_1[1]!=mutation_2[1]:
        # can't match unless they are the same type
        return False
    else:
        # if they are of the same type
        mutation_type = mutation_1[1]
        # do different things based on what type of mutation
        if mutation_type == INDEL:
            if fabs(mutation_1[2]-mutation_2[2]) <= INDEL_EDGE_TOLERANCE: 
                return True
            else:
                return False
                
        elif mutation_type == REPEAT:
            key_1, dummy1, position_1, strand_1, repeat_name_1, repeat_end_1 =   mutation_1
            key_2, dummy2, position_2, strand_2, repeat_name_2, repeat_end_2 = mutation_2
    
            if (fabs(position_1-position_2) <= REPEAT_EDGE_TOLERANCE) and (strand_1==strand_2) and (repeat_name_1==repeat_name_2):
                if version!='long':
                    # don't count orientation of IS element as unique
                    return True
                else:
                    if repeat_end_1==repeat_end_2:
                        return True
                    else:
                        return False
                        
        else:
            key_1, dummy1, position_L1, strand_L1, position_R1, strand_R1 = mutation_1
            key_2, dummy2, position_L2, strand_L2, position_R2, strand_R2 = mutation_2

            if fabs(position_L1-position_L2) <= OTHER_EDGE_TOLERANCE and strand_L1==strand_L2 and fabs(position_R1-position_R2) <= OTHER_EDGE_TOLERANCE and strand_R1==strand_R2:
                return True
            else:
                return False
                    
def get_registered_mutation(mutation_string,version='long'):
    return str_mutation_map[version][mutation_string]

def get_registered_mutation_string(mutation, version='long'):
    
    for registered_mutation in mutation_str_map[version].keys():
        if fuzzy_match(mutation, registered_mutation,version):
            #sys.stderr.write("%s %s\n" % (str(mutation), str(registered_mutation)))
            return mutation_str_map[version][registered_mutation]

    # no fuzzy match,
    # so create new entry

    mutation_type = mutation[1]

    if mutation_type == INDEL:
        long_mutation_string = "indel_%d" % mutation[2]
        short_mutation_string = long_mutation_string
        indel_alleles[short_mutation_string] = set()

    elif mutation_type == REPEAT:
        long_mutation_string = "junction_%d_%d_%s_%s" % (mutation[2],mutation[3],mutation[4],mutation[5])
        short_mutation_string = "junction_%d_%d_%s" % (mutation[2], mutation[3], mutation[4]) 
    else:
        long_mutation_string = "junction_%d_%d_%d_%d" % (mutation[2],mutation[3],mutation[4],mutation[5])
        short_mutation_string = long_mutation_string
    
    if version=='long':
        mutation_string = long_mutation_string
    else:
        mutation_string = short_mutation_string

    mutation_str_map[version][mutation] = mutation_string
    str_mutation_map[version][mutation_string] = mutation
    return mutation_string

# parse junctions from rebreseq files
for t,sample_name in zip(sample_times, sample_list):
    #print t, sample_name, gd_filename_map[sample_name]
    if not t in mutation_map:
        mutation_map[t] = {}
    gd_file = open(gd_filename_map[sample_name],"r")
    sys.stderr.write("%s\n" % sample_name)

    initial_sv_list = []
    initial_indel_list = []

    for line in gd_file:

        # we are only looking at JC evidence here
        if not line.startswith('JC'):
            continue
 
        items = line.split()
        unnamed_items = []
        named_items = {}
        for item in items:
            if '=' in item:
                subitems = item.split('=')
                named_items[subitems[0]] = subitems[1]
            else:
                unnamed_items.append(item)

        key = named_items['key']

        if named_items['side_1_annotate_key'] == 'repeat' and named_items['side_2_annotate_key'] == 'repeat':
            # we can't do anything about junctions with both ends in repeat regions
            continue
        
        elif named_items['side_1_annotate_key'] == 'repeat':
            
            # side 1 is a repeat and side 2 is not
            involves_repeat = True
            # get gene position from side 2
            position = long(unnamed_items[7]) 
            # get strand position from side 2
            strand = long(unnamed_items[8])
            R = float(named_items['side_2_read_count'])
            wR = 1.0/(1+float(named_items['side_2_possible_overlap_registers']))

            # repeat stuff 
            repeat_position = long(unnamed_items[4])
            idx = get_closest_repeat_idx(repeat_position, repeat_data)

            repeat_name = repeat_names[idx].split(":")[1]
            repeat_complement = repeat_complements[idx]

            # find out whether we are near start or finish
            distances = numpy.fabs( numpy.array([repeat_starts[idx], repeat_ends[idx]]) - repeat_position )
            is_end = distances.argmin()
            if repeat_complement:
                is_end = not is_end
            if is_end:
                repeat_end = 'end'
            else:
                repeat_end = 'start'
 
            mutation = (key, REPEAT, position, strand, repeat_name, repeat_end)
        
            A = float(named_items['new_junction_read_count'])
            wA = 1.0/(1+float(named_items['junction_possible_overlap_registers']))

            initial_sv_list.append( (mutation, wA, A, wR, R) )
        
        
        elif named_items['side_2_annotate_key'] == 'repeat':
        
            # side 2 is a repeat and side 1 is not
        
            # get position from side 1
            position = long(unnamed_items[4])
            # get strand from side 1
            strand = long(unnamed_items[5])
            R = float(named_items['side_1_read_count'])
            wR = 1.0/(1+float(named_items['side_1_possible_overlap_registers']))

            # repeat stuff 
            repeat_position = long(unnamed_items[7])
            idx = get_closest_repeat_idx(repeat_position, repeat_data)

            repeat_name = repeat_names[idx].split(":")[1]
            repeat_complement = repeat_complements[idx]

            # find out whether we are near start or finish
            distances = numpy.fabs( numpy.array([repeat_starts[idx], repeat_ends[idx]]) - repeat_position )
            is_end = distances.argmin()
            if repeat_complement:
                is_end = not is_end
            if is_end:
                repeat_end = 'end'
            else:
                repeat_end = 'start'

            mutation = (key, REPEAT, position, strand, repeat_name, repeat_end)
            
            A = float(named_items['new_junction_read_count'])
            wA = 1.0/(1+float(named_items['junction_possible_overlap_registers']))

            initial_sv_list.append( (mutation, wA, A, wR, R) )
            
        else:
            # neither side is a repeat
            position_1 = long(unnamed_items[4])
            strand_1 = long(unnamed_items[5])
            position_2 = long(unnamed_items[7])
            strand_2 = long(unnamed_items[8])
            R1 = float(named_items['side_1_read_count'])
            w1 = 1.0/(1+float(named_items['side_1_possible_overlap_registers']))
            R2 = float(named_items['side_2_read_count'])
            w2 = 1.0/(1+float(named_items['side_2_possible_overlap_registers']))

            wR = w1+w2
            R = (w1*R1+w2*R2)/wR
            
            # make sure smaller position is position 1
            # (should already be true?)
            if position_2 < position_1:
                sys.stderr.write("Switching positions: %d:%d %d:%d\n" % (position_1, strand_1, position_2, strand_2))

                position_1, position_2 = position_2, position_1
                strand_1, strand_2 = strand_2, strand_1
            
            A = float(named_items['new_junction_read_count'])
            wA = 1.0/(1+float(named_items['junction_possible_overlap_registers']))

            if strand_1==-1 and strand_2==1 and (position_2-position_1) <= INDEL_LENGTH_TOLERANCE:
                
                # simple indel
                if position_2 - position_1 > 1:
                    # deletion
                    mutation = (key, INDEL, position_1, '-%d' % (position_2-position_1-1))
                else:
                    # insertion
                    mutation = (key, INDEL, position_1, '+%s' % named_items['unique_read_sequence'])
               
                initial_indel_list.append((mutation, A, R))
            elif strand_1==1 and strand_2==-1 and (position_2-position_1) <= INDEL_LENGTH_TOLERANCE:
                
                # a duplication

                duplication_length = position_2-position_1
                if 'unique_read_sequence' in named_items:
                    duplication_length += len(named_items['unique_read_sequence'])                

                mutation = (key, INDEL, position_1-1, '+%d' % duplication_length)
                initial_indel_list.append((mutation,A,R))

            else:     
                # not a simple indel, most likely inversion
                mutation = (key, OTHER, position_1, strand_1, position_2, strand_2)
                initial_sv_list.append( (mutation, wA, A, wR, R) )

    gd_file.close()


    # condense indel list
    merged_indel_map = {}
    for mutation, A, R in initial_indel_list:

        # don't bias averages from things that 
        # truly are not there        
        # Q: for indels, how do these things even get here?
        if A+R < 1:
            continue
 
        mutation_string = get_registered_mutation_string(mutation)    
        #if mutation_string == 'indel_4558685':
        #    sys.stderr.write("%s %s\n" % (mutation_string, str(mutation)))
        
        if mutation_string not in merged_indel_map:
            merged_indel_map[mutation_string] = (0,0,0)
        
        indel_alleles[mutation_string].add(mutation[3])

        old_n,old_A, old_R = merged_indel_map[mutation_string]
        new_n = old_n+1
        new_A = old_A+A 
        # we have a choice about how to do this
        # technically they should all be the same,
        # but in practice they differ for a few reasons
        # could either average or max
        new_R = max([old_R,R])
        #new_R = (old_R*old_n+R)/new_n
        
        merged_indel_map[mutation_string] = (new_n,new_A, new_R)
            
    # now condense sv list
    merged_mutation_map = {}
    for mutation, wA, A, wR, R in initial_sv_list:

        # if no reads, don't even include it! 
        if A+R < 1:
            continue
        
        mutation_string = get_registered_mutation_string(mutation)
        if mutation_string not in merged_mutation_map:
            merged_mutation_map[mutation_string] = (0,0,0,0,0)
       
        # essentially averages things together with 1/registers as weights 
        old_n, old_wA, old_A, old_wR, old_R = merged_mutation_map[mutation_string]
        new_n = old_n+1
        new_wA = old_wA + wA
        new_A = (old_wA*old_A/(old_n+(old_n==0)) + wA*A) / (new_wA) * new_n
        new_wR = old_wR + wR
        new_R = (old_wR*old_R + wR*R) / (new_wR) 
        merged_mutation_map[mutation_string] = (new_n,new_wA, new_A, new_wR, new_R)


    # account for weird normalization issues for non-indels
    renormalized_mutation_map = {}
    for mutation_string in merged_mutation_map.keys():
        
        
        short_mutation_string = get_registered_mutation_string(get_registered_mutation(mutation_string),version='short')
        
        #sys.stderr.write("%s\n" % short_mutation_string)
        
        if short_mutation_string not in renormalized_mutation_map:
            renormalized_mutation_map[short_mutation_string] = {}
        renormalized_mutation_map[short_mutation_string][mutation_string] = merged_mutation_map[mutation_string]
        
        #sys.stderr.write("%s %s\n" % (mutation_string, short_mutation_string))

    for short_mutation_string in renormalized_mutation_map.keys(): 
        # calculate renormalized wD, D
        
        # first count the total n,A, and wA
        # across all alternate alleles with the same short mutation string
        # i.e., the same except for switching the orientation of the IS element
        wAtot = 0
        Atot = 0
        ntot = 0 
        for mutation_string in renormalized_mutation_map[short_mutation_string].keys():
            n,wA,A,wR,R = renormalized_mutation_map[short_mutation_string][mutation_string]
            ntot += n
            wAtot += wA
            Atot += A*wA/n 
            #if len(renormalized_mutation_map[short_mutation_string].keys()) > 1:
            #    sys.stderr.write("%g %g\n" % (A,R))
        
        # this is the total number of alternate alleles
        Atot = Atot/wAtot*ntot

        #if len(renormalized_mutation_map[short_mutation_string].keys()) > 1:
        #    sys.stderr.write("Done!\n")

        #if len(renormalized_mutation_map[short_mutation_string].keys()) > 1:
            #sys.stderr.write("%s\n" % (" ".join(renormalized_mutation_map[short_mutation_string].keys())))

        for mutation_string in renormalized_mutation_map[short_mutation_string].keys():
            n,wA,A,wR,R = renormalized_mutation_map[short_mutation_string][mutation_string]
            location = long(mutation_string.split("_")[1])
            alt = A

            #       total alt weight           
            #                  |
            depth = (Atot*wAtot/ntot + R*wR/n)/(wA/n)
            #sys.stderr.write("%s %g %g\n" % (mutation_string, alt, depth)) 

            if mutation_string in mutation_map[t]:
                # independently sequenced same timepoint: add coverage up
                old_mutation, old_alt, old_depth = mutation_map[t][mutation_string]
                alt = old_alt + alt
                depth = old_depth + depth    

            mutation_map[t][mutation_string] = (("REL606", location, mutation_string),  alt, depth)
            
    for indel_string in merged_indel_map.keys(): 
        location = long(indel_string.split("_")[1])
        n,A,R = merged_indel_map[indel_string]
        D = A+R
        if indel_string in mutation_map[t]:
            # independently sequenced same timepoint: add coverage up
            old_mutation, old_A, old_D = mutation_map[t][indel_string]
            A = old_A + A
            D = old_D + D
        mutation_map[t][indel_string] = (("REL606", location, indel_string),  A, D)

            
mutation_timecourse_map = {}
for t in mutation_map.keys():
    for mutation,alt,depth in mutation_map[t].values():
        if mutation not in mutation_timecourse_map:
            mutation_timecourse_map[mutation] = {}
        mutation_timecourse_map[mutation][t] = (alt,depth)

rebreseq_junctions = {}
rebreseq_indels = {}

for mutation in sorted(mutation_timecourse_map.keys(), key=itemgetter(1)):
    times = []
    alts = []
    depths = []
    for t in sample_times_set:
        times.append(t)
        if t in mutation_timecourse_map[mutation]:
            alt = mutation_timecourse_map[mutation][t][0]
            depth = mutation_timecourse_map[mutation][t][1]
        else:
            alt = 0
            depth = 0
        alts.append(alt)
        depths.append(depth)

    chromosome,location,allele_str = mutation 
    times = numpy.array(times)
    alts = numpy.array(alts)
    depths = numpy.array(depths)
    if allele_str.startswith("junction"):
        rebreseq_junctions[(chromosome,location,allele_str)] = (chromosome,location,allele_str,times,alts,depths)
    else:
        # a simple indel
        rebreseq_indels[(chromosome,location)] = (chromosome,location,"indel;"+(";".join(subitem for subitem in indel_alleles[allele_str])),times,alts,depths)

#if ('REL606', 4558685) in rebreseq_indels:
#    sys.stderr.write("In rebreseq indels!\n")
#else:
#    sys.stderr.write("Not in rebreseq indels!\n")

# parse through trajectories from mpileup file
for line in sys.stdin:
    items = line.split(",")
    chromosome = items[0].strip()
    location = long(items[1])
    allele = items[2].strip()
    if not allele.startswith('indel'):
        # a snp, let it pass
        print line,
    else:
        if not (chromosome,location) in rebreseq_indels:
            # an indel without any additional things to add
            # let it pass
            print line,    
        else:
            # it's in the rebreseq list!
            
            # parse and add the things to it! 
            times = numpy.array([float(subitem) for subitem in items[3].split()])
            alts = numpy.array([float(subitem) for subitem in items[4].split()])
            depths = numpy.array([float(subitem) for subitem in items[5].split()])

            other_allele = rebreseq_indels[(chromosome,location)][2]
            other_alts = rebreseq_indels[(chromosome,location)][4]
            other_depths = numpy.ceil(rebreseq_indels[(chromosome,location)][5])

            # debug stuff
            #if location==4375766:
            #    for time_idx in xrange(0,len(times)):
            #        sys.stderr.write("%d: %d\n" % (location, times[time_idx]))
            #        sys.stderr.write("Am = %g, Dm = %g\n" % (alts[time_idx], depths[time_idx]))
            #        sys.stderr.write("Ab = %g, Db = %g\n" % (other_alts[time_idx], other_depths[time_idx]))
            


            old_alleles = allele.split(";")
            new_alleles = other_allele.split(";")
           
            total_alleles = set(old_alleles[1:])
            total_alleles.update(new_alleles[1:])
            
            new_allele = "indel;"+(";".join(total_alleles))

            new_times = times
            new_alts = alts+other_alts
            # Two ways to do this:
            # add junction alts to mpileup depths:
            # new_depths = depths + other_alts
            # this typically double counts junction alts, but in weird ways

            # or mpileup alts to breseq depths:
            #new_depths = other_depths + alts
            # can't tell, but maybe double counts breseq alts in weird ways?

            # or just take breseq depths
            new_depths = other_depths*(other_alts>0.5)+depths*(other_alts<0.5)
 
            # however, this doesn't always work, so print out the violators 
            if (new_alts > new_depths).any():
                sys.stderr.write("New alts exceed depths: %d\n" % location)
                sys.stderr.write("%s\n" % new_allele)
            
                for i in xrange(0,len(times)):
                    if new_alts[i]-new_depths[i] > 1:
                        sys.stderr.write("%g: (%g,%g) (%g,%g)\n" % (times[i],alts[i],depths[i],other_alts[i],other_depths[i]))                   

                new_depths = numpy.fmax(new_depths,new_alts)

            # print it
            time_strs = ["%g" % t for t in new_times]
            alt_strs = ['%g' % a for a in new_alts]
            depth_strs = ['%0.1f' % d for d in new_depths]
            print ", ".join([chromosome, str(location), new_allele, " ".join(time_strs), " ".join(alt_strs), " ".join(depth_strs)])

            # delete it from map
            del rebreseq_indels[(chromosome,location)]

# done parsing through that!

output_list = []

# now print the rest
for chromosome,location,allele,times,alts,depths in rebreseq_indels.values():    
    if (alts >= 2).sum() > 1 and ((alts>=2)*(depths>=10)*(alts>=(0.05*depths))).sum()>0: # recording condition
        time_strs = ["%g" % t for t in times]
        alt_strs = ['%g' % a for a in alts]
        depth_strs = ['%0.1f' % d for d in depths]
        output_string = ", ".join([chromosome, str(location), allele, " ".join(time_strs), " ".join(alt_strs), " ".join(depth_strs)])
        output_list.append((chromosome,location,output_string))

# now print the rest
for chromosome,location,allele,times,alts,depths in rebreseq_junctions.values():
    if (alts >= 2).sum() > 1 and ((alts>=2)*(depths>=10)*(alts>=(0.05*depths))).sum()>0: # recording condition
        time_strs = ["%g" % t for t in times]
        alt_strs = ['%g' % a for a in alts]
        depth_strs = ['%0.1f' % d for d in depths]
        output_string =  ", ".join([chromosome, str(location), allele, " ".join(time_strs), " ".join(alt_strs), " ".join(depth_strs)])
        output_list.append((chromosome,location,output_string))

output_list.sort()
for chromosome,location,output_string in output_list:
    print output_string
