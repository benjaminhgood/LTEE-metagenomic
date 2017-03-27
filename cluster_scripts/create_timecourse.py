# This script takes the output from an mpileup of our 
# second breseq run (the "rebreseq" step") and calculates 
# a timecourse for mutations at each site in the genome. 
#
# Distinct SNPs at the same site are treated as two separate mutations
# but indels are merged together into a single alternate allele

import numpy
import sys
import population_parameters

reference_strs = set(['.', ','])
snp_strs = set(['A','C','T','G', 'N','*'])
indel_strs = set(['+','-'])
refskip_strs = ['<','>']

mpileup_file = sys.stdin

first_position = long(sys.argv[1])
last_position = long(sys.argv[2])
population = sys.argv[3]
sample_names = sys.argv[4:]

sample_list = []
sample_times = []

for sample_name, sample_time in zip(population_parameters.sample_names[population], population_parameters.sample_times[population]):
    if sample_name in sample_names:
        sample_list.append(sample_name)
        sample_times.append(sample_time)

num_lines = 0
num_printed_lines = 0

for line in mpileup_file:
    num_lines += 1
    if num_lines % 1000 == 0:
        sys.stderr.write("%d lines processed, %d passed\n" % (num_lines,num_printed_lines))
        #sys.exit(0)

    items = line.split('\t')
    chromosome = items[0]
    position = long(items[1])
    ref = items[2]
 
    if position < first_position:
        continue
    
    if position > last_position:
        break

    times = []
    alts = []
    depths = []
    indel_depth_reductions = []

    num_samples = (len(items)-3)/3

    alleles = {ref: [0]*num_samples} 
    
    indel_alleles = set()
    
    for i in xrange(0, num_samples):
        t = sample_times[i]
            
        depth = long(items[(i+1)*3])
        allele_str = items[(i+1)*3+1].upper()
        depths.append(depth)
        indel_depth_reductions.append(0)
        times.append(t)
        j=0
        jmax = len(allele_str)

        # walk through the mpileup string and parse the alleles
        while j < jmax:
            if allele_str[j] == '^': # beginning of read
                j+=2 # next character is map quality, so skip it too

            elif j+1<jmax and (allele_str[j+1] in indel_strs):
                
                #sys.stderr.write("%s\n" % allele_str[j:])
 
                # an indel

                # only record the base allele if different from ref
                if allele_str[j] in reference_strs:
                    base_allele=''
                else:
                    base_allele=allele_str[j]

                j+=1

                # whether it is a plus or minus
                indel_allele = allele_str[j]
                j+=1

                # calculate length of insertion or deletion
                # (i.e., how many characters to check afterwards)
                k = j
                while allele_str[k].isdigit():
                    k+=1
                indel_len = long(allele_str[j:k])
                j = k

                # the inserted or deleted bases themselves
                indel_bases = allele_str[j:j+indel_len]
                j += indel_len

                if indel_allele=='+':
                    # insertion
                    full_indel_allele = ('%s+%s' % (base_allele, indel_bases))
                else:
                    # deletion
                    full_indel_allele = ('%s-%d' % (base_allele, indel_len))                
                indel_alleles.add(full_indel_allele)
                
                if 'indel' not in alleles:
                    alleles['indel'] = [0]*num_samples
                alleles['indel'][i] += 1

            elif allele_str[j] in reference_strs:
                # reference

                alleles[ref][i] += 1
                j+=1

                if j<jmax and allele_str[j]=='$':
                    # ref fell at end of read
                    # don't count it for support for indel
                    indel_depth_reductions[i] += 1
                    j+=1
               
            elif allele_str[j] in snp_strs: # regular SNP
                
                if allele_str[j] not in alleles:
                    alleles[allele_str[j]] = [0]*num_samples
                alleles[allele_str[j]][i] += 1
                j+=1

            else:
                # not sure, do nothing
                j+=1 
    
    depths = numpy.array(depths) 
    indel_depth_reductions = numpy.array(indel_depth_reductions)

    depth_map = {}
    indel_depth_map = {}
    for i in xrange(0,len(depths)):
        if times[i] not in depth_map:
            depth_map[times[i]]=0
            indel_depth_map[times[i]]=0
        depth_map[times[i]] += depths[i]
        indel_depth_map[times[i]] += (depths[i]-indel_depth_reductions[i])

    merged_times = numpy.array([t for t in sorted(depth_map.keys())])
    merged_depths = numpy.array([depth_map[t] for t in merged_times])
    merged_indel_depths = numpy.array([indel_depth_map[t] for t in merged_times])
    alt_map = {}

    for allele in alleles.keys():
         
        # don't do anything about refs
        if allele==ref or allele=='*':
            continue
       
        allele_key = allele

        if allele_key not in alt_map:
            alt_map[allele_key] = {t: 0 for t in depth_map.keys()}

        for i in xrange(0,len(depths)):
            alt_map[allele_key][times[i]] += alleles[allele][i]

    # get merged alts for each allele (max 5)
    merged_alts = {}
    for allele_key in alt_map.keys():
        merged_alts[allele_key] = numpy.array([alt_map[allele_key][t] for t in merged_times])
    
    if 'indel' in merged_alts.keys():
        # if basically an indel
        if ((merged_alts['indel']*1.0/(merged_indel_depths+(merged_indel_depths==0)) > 0.5)*(merged_indel_depths>10)).any():
            # merge everything together
            for allele_key in merged_alts.keys():
                if allele_key=='indel':
                    continue
                
                indel_alleles.add(allele_key)
                merged_alts['indel'] += merged_alts[allele_key]
            # delete reference to other SNPs
            merged_alts = {'indel' : merged_alts['indel']}
                        

    for allele_key in merged_alts.keys():
        if allele_key == 'indel':
            alt_allele = "indel;" + ";".join(indel_alleles)
            allele_depths = merged_indel_depths
        else:
            alt_allele = "%s->%s" % (ref,allele_key)
            allele_depths = merged_depths

        allele_alts = merged_alts[allele_key]

        if ((allele_alts>=2).sum() > 1) and ((allele_alts >= 2)*(allele_depths>=10)*(allele_alts >= 0.05*allele_depths)).sum() > 0:
            print ", ".join([chromosome, str(position), alt_allele, " ".join(str(t) for t in merged_times), " ".join(str(a) for a in allele_alts), " ".join(str(d) for d in allele_depths)])
            num_printed_lines+=1
            
            
