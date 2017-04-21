import timecourse_utils
import parse_file
from math import log,exp
import numpy
from scipy.special import gammaln as loggamma
import sys
from scipy.optimize import fmin   
from scipy.interpolate import interp1d
from scipy.stats import nbinom     

import uniform_prior_calculate_haplotypes as clade_model

if len(sys.argv)>1:
    populations = sys.argv[1:]
else:
    populations = parse_file.all_lines

min_coverage = 5
#min_clade_length = 0
min_clade_length = 9750

states = parse_file.clade_hmm_states

for population in populations:
 
    mutations = []
    good_mutations = []
 
    sys.stderr.write("Processing %s\n" % population)
 
    mutation_data, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    times = mutation_data[0][10]
    
    for mutation_idx in xrange(0,len(mutation_data)):
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutation_data[mutation_idx] 
        
        good_idxs, masked_alts, masked_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue, min_coverage)
 
        max_freq = (masked_alts*1.0/(masked_depths+(masked_depths<1))).max()
        num_good_timepoints = (masked_depths>0.5).sum()
        
        if num_good_timepoints > 5:
            mutations.append((masked_alts*1.0, masked_depths*1.0))
            if var_type!='sv' and var_type!='indel' and masked_depths[-1]>1:
                good_mutations.append((masked_alts*1.0, masked_depths*1.0))
        else:
            mutations.append((numpy.zeros_like(times)*1.0, numpy.zeros_like(times)))
    
    # run first round with just good mutations
    # i.e., SNPs that last until the end
    A = []
    D = []
    for i in xrange(0,len(good_mutations)):
        Ai,Di = good_mutations[i]
        A.append(Ai)
        D.append(Di)
        
    A = numpy.array(A)
    D = numpy.array(D)
    
    # Run HMM decode
    # estimates clade frequencies from scratch
    Pstate, Ls, fmajors, fminors = clade_model.infer_hmm(A,D,num_iterations=10)
    #Pstate, fmajors, fminors = clade_model.infer_hmm(A,D,num_iterations=1)
    #Pstate, fmajors, fminors = clade_model.infer_hmm(A,D,num_iterations=1,infer_fs=False)
    
    
    # calculate "length" of coexistance 
    nonzero_idxs = numpy.nonzero(((fmajors+fminors)>0.8)*((2*fmajors*(1-fmajors))>0.32)*((2*fminors*(1-fminors))>0.32))[0]
    
    if len(nonzero_idxs)==0:
        length=0
    else:
        endpoint = times[nonzero_idxs[-1]]
        startpoint = times[numpy.nonzero((fmajors+fminors)>0.8)[0][0]]
        length = endpoint-startpoint
    
    
    if length>min_clade_length:
        # real frequency dependence!
        pass
    else:
        # not real frequency depencence 
        # e.g., maybe just clonal interference
        # set the clades to zero
        fminors = numpy.zeros_like(fminors)
        fmajors = numpy.zeros_like(fmajors)
    
    sys.stdout.write("%s: length of coexistence = %g\n" % (population,length))

        
    # now run again with all mutations
    # (but treat clade frequencies as fixed)
    A = []
    D = []
    for i in xrange(0,len(mutations)):
        Ai,Di = mutations[i]
        A.append(Ai)
        D.append(Di)   
    A = numpy.array(A)
    D = numpy.array(D)
    
    Pstate , Ls, fmajors, fminors = clade_model.infer_hmm(A,D,f0s=(fmajors,fminors),num_iterations=5,infer_fs=False)              
    
    # make thing that says if filtered:
    
    filtered_mutations = numpy.ones_like(A)
    for i in xrange(0,A.shape[0]):
        if D[i,:].sum() < 0.5:
            filtered_mutations[i,:] = numpy.zeros_like(times)
    
    ns = numpy.zeros((Pstate.shape[2],Pstate.shape[1]))
    
    for l in xrange(0,ns.shape[0]):
        ns[l,:] = (Pstate[:,:,l]*filtered_mutations).sum(axis=0)
    
    
    #Ls = numpy.zeros_like(A) 
    #for i in xrange(0,A.shape[0]):
    #    if filtered_mutations[i,0] < 0.5:
    #        Ls[i] = numpy.zeros_like(times)
    #    else:
    #        Ls[i] = (Pstate[i,:,:]).argmax(axis=1)
    
    hard_ns = numpy.zeros_like(ns) 
    for l in xrange(0,ns.shape[0]):
        hard_ns[l,:] = (Ls==l).sum(axis=0)
    
    #fmajors = fs*((hard_ns[5,:]+hard_ns[4,:])>0.5)
    #fminors = (1-fs)*((hard_ns[5,:]+hard_ns[4,:])>0.5)
           
    # Write stuff
    haplotype_filename = parse_file.data_directory+("%s_haplotype_timecourse.txt" % population)
    haplotype_file = open(haplotype_filename, "w")
    haplotype_file.write("%s\n" % ", ".join([str(t) for t in times]))
    haplotype_file.write("%s\n" % ", ".join([str(f) for f in fmajors]))
    haplotype_file.write("%s\n" % ", ".join([str(f) for f in fminors]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[0,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[5:,:].sum(axis=0)]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[1,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[2,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[3,:]]))
    haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[4,:]])) 
    for i in xrange(0,Ls.shape[0]):  
    #    start_idx=-25
    #    end_idx=-15
    #    if mutation_data[i][0]==4567083:
    #        print ", ".join(["%0.1f" % (t/1000.0) for t in times[start_idx:end_idx]])
    #        for state in parse_file.clade_hmm_states.keys():
    #            print state,  ", ".join(["%0.1e" % p for p in Pstate[i,start_idx:end_idx,parse_file.clade_hmm_states[state]]])
    #        print ", ".join(["%d" % a for a in A[i,start_idx:end_idx]])
    #        print ", ".join(["%d" % d for d in D[i,start_idx:end_idx]])
    #        print "Obs", ", ".join(["%0.2f" % f for f in (A*1.0/(D+(D==0)))[i,start_idx:end_idx]])
    #        print "Maj", ", ".join(["%0.2f" % f for f in fmajors[start_idx:end_idx]]) 
    #        print "Min", ", ".join(["%0.2f" % f for f in fminors[start_idx:end_idx]])   
        haplotype_file.write("%s\n" % ", ".join([str(l) for l in Ls[i,:]]))
        
    haplotype_file.close()
    