import sys
import pylab
import numpy
import bz2
import parse_file

population = sys.argv[1]
desired_locations = set([long(item) for item in sys.argv[2:]])

pylab.figure(1,figsize=(17,2))
pylab.figure(2,figsize=(17,2))
pylab.figure(3,figsize=(17,2))


mutations, depth_tuple = parse_file.parse_annotated_timecourse(population, only_passed=False)

population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple  
    
    

for mutation_idx in xrange(0,len(mutations)):
 
    location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
    if location in desired_locations:
            
        depth_ratios = depths*1.0/population_avg_depths
        
        times = times[population_avg_depths>5]
        alts = alts[population_avg_depths>5]
        depth_ratios = depth_ratios[population_avg_depths>5]
        depths = depths[population_avg_depths>5]
        
        freqs = (alts*1.0/(depths+(depths==0)))
        
        clone_depth_ratios = clone_depths*1.0/clone_avg_depths
        clone_times = clone_times[clone_avg_depths>5]
        clone_alts = clone_alts[clone_avg_depths>5]
        clone_depth_ratios = clone_depth_ratios[clone_avg_depths>5]
        
        clone_freqs = (clone_alts*1.0/(clone_depths+(clone_depths==0)))
        
        print clone_freqs
        
        pylab.figure(1)
        line, = pylab.plot(times,freqs,'.-')
        pylab.plot(clone_times,clone_freqs,'s',color=pylab.getp(line,'color'))
        
        pylab.figure(2)
        pylab.semilogy(times,depths+0.1*(depths==0),'.-',color=pylab.getp(line,'color'))
        pylab.semilogy(clone_times,clone_depths+0.1*(clone_depths==0),'s',color=pylab.getp(line,'color'))
        pylab.figure(3)
        pylab.semilogy(times,depth_ratios,'.-',color=pylab.getp(line,'color'))
        pylab.semilogy(clone_times,clone_depth_ratios,'s',color=pylab.getp(line,'color'))
        pylab.figure(4)
        pylab.semilogx(depth_ratios, freqs,'.',color=pylab.getp(line,'color'))
        pylab.semilogx(clone_depth_ratios, clone_freqs,'s',color=pylab.getp(line,'color'))
        
        
pylab.figure(1)
#pylab.xlim([0,64])
pylab.ylim([-0.05,1.05])
pylab.figure(2)
#pylab.xlim([0,64])
pylab.plot([0,60000],[5,5],'k:')
pylab.figure(3)
#pylab.xlim([0,64])
pylab.figure(4)
pylab.ylim([-0.1,1.1])
print "Done!"
pylab.show()