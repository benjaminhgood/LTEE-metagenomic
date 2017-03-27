import sys
import numpy
from math import log10,ceil
import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import parse_file
import timecourse_utils
from scipy.stats import fisher_exact
import statsmodels.api as sm
import mutation_spectrum_utils
from scipy.stats import beta
from numpy.random import shuffle
import figure_utils

theory_times = numpy.arange(0,125)*500

metapopulations = ['nonmutator','mutator']
metapopulation_labels = {'nonmutator': 'Nonmutators', 'mutator':'Mutators'}
metapopulation_populations = {'nonmutator': parse_file.complete_nonmutator_lines, 'mutator': parse_file.mutator_lines}

colors = {'nonmutator': figure_utils.nonmutator_group_color, 'mutator': figure_utils.mutator_group_color}


##############################################################################
#
# Set up figures
#
##############################################################################

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

pfix_axes = {}
parallelism_axes = {}

##############################################################################
#
# First set up Fig. 5 (rates of fitness and mutation accumulation with time)
#
##############################################################################

fig5 = plt.figure(figsize=(7.2, 4.5))

outer_grid  = gridspec.GridSpec(2, 1, height_ratios=[2,0.9], hspace=0.20)

inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, width_ratios=[3.4,1], height_ratios=[1,1],hspace=0.10,
                subplot_spec=outer_grid[0], wspace=0.20) #, hspace=0.08)

pfix_grid = gridspec.GridSpecFromSubplotSpec(1, 3,
                width_ratios=[1.3,0.9,0.9],
                wspace=0.10,
                subplot_spec=outer_grid[1])


###########
#
# Nonmutator variant types through time
#
###########
nonmutator_axis = plt.Subplot(fig5, inner_grid[0,0])
fig5.add_subplot(nonmutator_axis)

nonmutator_axis.spines['top'].set_visible(False)
nonmutator_axis.spines['right'].set_visible(False)
nonmutator_axis.get_xaxis().tick_bottom()
nonmutator_axis.get_yaxis().tick_left()
nonmutator_axis.get_xaxis().set_tick_params(direction='out')
nonmutator_axis.set_title('Nonmutators', y=0.85,fontsize=5)

###########
#
# Mutator variant types through time
#
###########
mutator_axis = plt.Subplot(fig5, inner_grid[1,0])
fig5.add_subplot(mutator_axis)
mutator_axis.set_ylabel('Total mutations detected ($\\times 10^3$)',y=1.1)
mutator_axis.set_xlabel('Generation, $t$')

mutator_axis.spines['top'].set_visible(False)
mutator_axis.spines['right'].set_visible(False)
mutator_axis.get_xaxis().tick_bottom()
mutator_axis.get_yaxis().tick_left()
mutator_axis.set_title('Mutators',y=0.85,fontsize=5)

#mutator_axis.get_yaxis().set_tick_params(direction='out')
#mutator_axis.get_xaxis().set_tick_params(direction='out')

metapopulation_axes = {'nonmutator':nonmutator_axis, 'mutator':mutator_axis}

###########
#
# Nonmutator multiplicity distribution
#
###########

nonmutator_parallelism_axis = plt.Subplot(fig5, inner_grid[0,1])

fig5.add_subplot(nonmutator_parallelism_axis)

nonmutator_parallelism_axis.spines['top'].set_visible(False)
nonmutator_parallelism_axis.spines['right'].set_visible(False)
nonmutator_parallelism_axis.get_xaxis().tick_bottom()
nonmutator_parallelism_axis.get_yaxis().tick_left()

###########
#
# Mutator multiplicity distribution
#
###########

mutator_parallelism_axis = plt.Subplot(fig5, inner_grid[1,1])
fig5.add_subplot(mutator_parallelism_axis)

mutator_parallelism_axis.set_ylabel('Fraction mutations $\geq m$',y=1.1)
mutator_parallelism_axis.set_xlabel('Gene multiplicity, $m$')

mutator_parallelism_axis.spines['top'].set_visible(False)
mutator_parallelism_axis.spines['right'].set_visible(False)
mutator_parallelism_axis.get_xaxis().tick_bottom()
mutator_parallelism_axis.get_yaxis().tick_left()

parallelism_axes[True] = {'nonmutator': nonmutator_parallelism_axis, 'mutator': mutator_parallelism_axis}

###########
#
# Fixation probability vs variant type
#
###########

vartype_pfix_axis = plt.Subplot(fig5, pfix_grid[0])
fig5.add_subplot(vartype_pfix_axis)

vartype_pfix_axis.spines['top'].set_visible(False)
vartype_pfix_axis.spines['right'].set_visible(False)
vartype_pfix_axis.get_xaxis().tick_bottom()
vartype_pfix_axis.get_yaxis().tick_left()

vartype_pfix_axis.set_ylabel('$\mathrm{Pr}[\mathrm{fixed}|\mathrm{detected}$')
vartype_pfix_axis.set_ylim([0,1.05])
vartype_pfix_axis.set_xlim([-0.5,len(parse_file.var_types)-0.5])
vartype_pfix_axis.set_xticks(numpy.arange(0,len(parse_file.var_types)))
vartype_pfix_axis.set_xticklabels(parse_file.var_types,rotation='vertical')

###########
#
# Nonmutator Fixation probability vs multiplicity 
#
###########

nonmutator_pfix_axis = plt.Subplot(fig5, pfix_grid[1])
fig5.add_subplot(nonmutator_pfix_axis)

nonmutator_pfix_axis.spines['top'].set_visible(False)
nonmutator_pfix_axis.spines['right'].set_visible(False)
nonmutator_pfix_axis.get_xaxis().tick_bottom()
nonmutator_pfix_axis.get_yaxis().tick_left()
nonmutator_pfix_axis.set_xlabel('Gene multiplicity, $m$')

###########
#
# Mutator Fixation probability vs multiplicity 
#
###########

mutator_pfix_axis = plt.Subplot(fig5, pfix_grid[2])
fig5.add_subplot(mutator_pfix_axis)

mutator_pfix_axis.spines['top'].set_visible(False)
mutator_pfix_axis.spines['right'].set_visible(False)
mutator_pfix_axis.get_xaxis().tick_bottom()
mutator_pfix_axis.get_yaxis().tick_left()
mutator_pfix_axis.set_xlabel('Gene multiplicity, $m$')

pfix_axes[True] = {'nonmutator': nonmutator_pfix_axis, 'mutator': mutator_pfix_axis}

##############################################################################
#
# Set up Fig. S10 (multiplicity plots without svs)
#
##############################################################################

nosv_fig = plt.figure(figsize=(5, 3.5))

nosv_outer_grid  = gridspec.GridSpec(2, 2, wspace=0.1, hspace=0.05)

###########
#
# Nonmutator multiplicity distribution (nosv)
#
###########

nonmutator_parallelism_axis_nosv = plt.Subplot(nosv_fig, nosv_outer_grid[0,0])

nosv_fig.add_subplot(nonmutator_parallelism_axis_nosv)

nonmutator_parallelism_axis_nosv.spines['top'].set_visible(False)
nonmutator_parallelism_axis_nosv.spines['right'].set_visible(False)
nonmutator_parallelism_axis_nosv.get_xaxis().tick_bottom()
nonmutator_parallelism_axis_nosv.get_yaxis().tick_left()

nonmutator_parallelism_axis_nosv.set_ylabel('Total mutations $\geq m$')
nonmutator_parallelism_axis_nosv.set_xlim([1,1e02])

###########
#
# Mutator multiplicity distribution (nosv)
#
###########

mutator_parallelism_axis_nosv = plt.Subplot(nosv_fig, nosv_outer_grid[0,1])
nosv_fig.add_subplot(mutator_parallelism_axis_nosv)

mutator_parallelism_axis_nosv.spines['top'].set_visible(False)
mutator_parallelism_axis_nosv.spines['right'].set_visible(False)
mutator_parallelism_axis_nosv.get_xaxis().tick_bottom()
mutator_parallelism_axis_nosv.get_yaxis().tick_left()
mutator_parallelism_axis_nosv.set_xlim([1,1e02])

parallelism_axes[False] = {'nonmutator': nonmutator_parallelism_axis_nosv, 'mutator': mutator_parallelism_axis_nosv}

###########
#
# Nonmutator Fixation probability vs multiplicity (nosv)
#
###########

nonmutator_pfix_axis_nosv = plt.Subplot(nosv_fig, nosv_outer_grid[1,0])
nosv_fig.add_subplot(nonmutator_pfix_axis_nosv)

nonmutator_pfix_axis_nosv.spines['top'].set_visible(False)
nonmutator_pfix_axis_nosv.spines['right'].set_visible(False)
nonmutator_pfix_axis_nosv.get_xaxis().tick_bottom()
nonmutator_pfix_axis_nosv.get_yaxis().tick_left()
nonmutator_pfix_axis_nosv.set_ylabel('Conditional fixation probability')
nonmutator_pfix_axis_nosv.set_xlabel('Gene multiplicity, $m$')
nonmutator_pfix_axis_nosv.set_ylim([0,1.05])
nonmutator_pfix_axis_nosv.set_xlim([1,1e02])


###########
#
# Mutator Fixation probability vs multiplicity 
#
###########

mutator_pfix_axis_nosv = plt.Subplot(nosv_fig, nosv_outer_grid[1,1])
nosv_fig.add_subplot(mutator_pfix_axis_nosv)

mutator_pfix_axis_nosv.spines['top'].set_visible(False)
mutator_pfix_axis_nosv.spines['right'].set_visible(False)
mutator_pfix_axis_nosv.get_xaxis().tick_bottom()
mutator_pfix_axis_nosv.get_yaxis().tick_left()
mutator_pfix_axis_nosv.set_xlabel('Gene multiplicity, $m$')
mutator_pfix_axis_nosv.set_ylim([0,1.05])
mutator_pfix_axis_nosv.set_xlim([1,1e02])

pfix_axes[False] = {'nonmutator': nonmutator_pfix_axis_nosv, 'mutator': mutator_pfix_axis_nosv}


#############################
#
# Now do actual calculations and plot figures
#
#############################
   

Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes()

for metapopulation in metapopulations:

    axis = metapopulation_axes[metapopulation]
    color = colors[metapopulation]
    label = metapopulation_labels[metapopulation]
    populations = metapopulation_populations[metapopulation]
    
    if metapopulation == 'nonmutator':
        metapopulation_offset = -0.35
    else:
        metapopulation_offset = 0

    #############################
    #
    # Variant types through time
    #
    #############################

    var_type_trajectories = {var_type: numpy.zeros_like(theory_times)*1.0 for var_type in parse_file.var_types}
    var_type_totals = {var_type: 0 for var_type in parse_file.var_types}
    var_type_fixeds = {var_type: 0 for var_type in parse_file.var_types}

    for population in populations:
    
        sys.stderr.write("Processing %s...\t" % parse_file.get_pretty_name(population))

        # calculate mutation trajectories
        # load mutations
        mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
        population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
        dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
        state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    
        num_processed_mutations = 0
        for mutation_idx in xrange(0,len(mutations)):
 
            num_processed_mutations+=1 
            
            position, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
            Ls = haplotype_trajectories[mutation_idx]
            state_Ls = state_trajectories[mutation_idx]
        
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
            freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
            masked_state_Ls = state_Ls[good_idxs]
            masked_Ls = Ls[good_idxs]
            
            t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls, masked_Ls)
            #fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1], masked_freqs[-1])
            fixed_weight = timecourse_utils.calculate_clade_fixed_weight(masked_Ls[-1], masked_freqs[-1])
        
            var_type_trajectories[var_type][theory_times>=t] += 1.0
            var_type_totals[var_type] += 1.0
            var_type_fixeds[var_type] += fixed_weight
        
        sys.stderr.write("added %d mutations!\n" % num_processed_mutations)

    
    var_type_trajectories = [var_type_trajectories[var_type] for var_type in parse_file.var_types]
    var_type_totals = numpy.array([var_type_totals[var_type] for var_type in parse_file.var_types])*1.0
    var_type_fixeds = numpy.array([var_type_fixeds[var_type] for var_type in parse_file.var_types])*1.0
    
    survival_probabilities = var_type_fixeds/var_type_totals
    pfix_upper = numpy.array([beta.ppf(0.84, var_type_fixeds[i]+1, var_type_totals[i]-var_type_fixeds[i]+1, loc=0, scale=1) for i in xrange(0,len(var_type_totals))])
    pfix_lower = numpy.array([beta.ppf(0.16, var_type_fixeds[i]+1, var_type_totals[i]-var_type_fixeds[i]+1, loc=0, scale=1) for i in xrange(0,len(var_type_totals))])
        
    vartype_pfix_axis.bar(numpy.arange(0,len(parse_file.var_types))+metapopulation_offset, survival_probabilities+0.05,width=0.35, bottom=-0.05, linewidth=0, facecolor=color,alpha=0.5)

    for i in xrange(0,len(parse_file.var_types)):
        if i==0:
            vartype_pfix_axis.bar([-10],[1],label=label,linewidth=0, facecolor=color,alpha=0.5)            
        vartype_pfix_axis.plot([i+metapopulation_offset+0.15]*2, [pfix_lower[i],pfix_upper[i]],color=color,alpha=0.5,linewidth=0.5)
    
    
    syn_total = var_type_totals[0]
    syn_fixed = var_type_fixeds[0]
    syn_pfix = syn_fixed*1.0/syn_total
    
    non_total = var_type_totals[1]+var_type_totals[2]
    non_fixed = var_type_fixeds[1]+var_type_fixeds[2]
    non_pfix = non_fixed*1.0/non_total
    
    oddsratio, pvalue = fisher_exact([[non_total, non_fixed], [syn_total, syn_fixed]])
    
    dnds_total = (non_total/syn_total)/(Lnon/Lsyn)
    dnds_fixed = (non_fixed/syn_fixed)/(Lnon/Lsyn)
    
    sys.stdout.write("Synonymous: Target=%g, Total=%0.1f, Fixed=%0.1f, pfix=%0.2f\n" % (Lsyn, syn_total, syn_fixed, syn_pfix))
    sys.stdout.write("Nonsynonymous: Target=%g, Total=%0.1f, Fixed=%0.1f, pfix=%0.2f\n" % (Lnon, non_total, non_fixed, non_pfix))
    sys.stdout.write("Total dN/dS=%0.2f\n" % (dnds_total))
    sys.stdout.write("Fixed dN/dS=%0.2f\n" % (dnds_fixed))
    sys.stdout.write("Pfix pvalue=%g\n" % pvalue)
    
    top = numpy.zeros_like(theory_times)*1.0
    bottom = numpy.zeros_like(theory_times)*1.0

    bottom_fixed = 0
    top_fixed = 0
    fixed_scale_factor = sum(var_type_totals)*1.0/sum(var_type_fixeds)

    for i in xrange(0,len(parse_file.var_types)):

        bottom = top*1.0
        top += var_type_trajectories[i]
    
        bottom_fixed = top_fixed*1.0
        top_fixed+=var_type_fixeds[i]
    
        color = figure_utils.get_var_type_color(parse_file.var_types[i])
    
        axis.fill_between(theory_times[:121], bottom[:121]/1000.0, top[:121]/1000.0, color=color)
        axis.fill_between([theory_times[120],theory_times[120]+2000],[bottom[120]/1000.0,bottom[120]/1000.0],[top[120]/1000.0,top[120]/1000.0],color=color)
        axis.bar([-1],[1],facecolor=color,edgecolor='none',label=parse_file.var_types[i])
        axis.fill_between([66000,69000], [bottom_fixed*fixed_scale_factor/1000.0,bottom_fixed*fixed_scale_factor/1000.0],[top_fixed*fixed_scale_factor/1000.0,top_fixed*fixed_scale_factor/1000.0], color=color)

    axis.set_xlim([0,70000])


    ###################################################################
    #
    # Gene parallelism analysis
    #
    ###################################################################

    for include_svs in [True, False]:

        if include_svs:
            convergence_matrix_filename = parse_file.data_directory + "gene_convergence_matrix.txt"
        else:
            convergence_matrix_filename = parse_file.data_directory + "gene_convergence_matrix_nosvs.txt"
            
        convergence_matrix = parse_file.parse_convergence_matrix(convergence_matrix_filename)

        # Calculate median appearance time
        pooled_appearance_times = []
        for gene_name in convergence_matrix.keys():
            for population in metapopulation_populations[metapopulation]:
                for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                    pooled_appearance_times.append(t)
                    
        tstar = numpy.median(pooled_appearance_times)
                    
        sys.stdout.write("Median appearance time = %g\n" % tstar)

        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,allowed_populations=metapopulation_populations[metapopulation],Lmin=100)
        
        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        sys.stdout.write("Total parallelism = %g (p=%g)\n" % (G,pvalue))
        
        predictors = []
        early_predictors = []
        late_predictors = []
        responses = []
    
        gene_hits = []
        gene_predictors = []
        gene_pfixs = []
    
        Ls = []
    
    
        for gene_name in convergence_matrix.keys():
        
            
            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']
        
            Ls.append(convergence_matrix[gene_name]['length'])
            m = gene_parallelism_statistics[gene_name]['multiplicity']
                  
            n = 0
            nfixed = 0
                    
            for population in metapopulation_populations[metapopulation]:
                for t,L,Lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                    
                    #fixed_weight = timecourse_utils.calculate_fixed_weight(L, f)
                    fixed_weight = timecourse_utils.calculate_clade_fixed_weight(Lclade,f)
                    
                    predictors.append(m)
                    responses.append(fixed_weight)
                    
                    n+=1
                    nfixed+=fixed_weight
                    
                    if t<=tstar:
                        early_predictors.append(m)
                    else:
                        late_predictors.append(m)
        
            if n > 0.5:
                gene_hits.append(n)
                gene_predictors.append(m)
                gene_pfixs.append(nfixed*1.0/n)
                
        Ls = numpy.array(Ls)  
        ntot = len(predictors)      
        mavg = ntot*1.0/len(Ls)
        
        predictors, responses = (numpy.array(x) for x in zip(*sorted(zip(predictors, responses), key=lambda pair: (pair[0])))) 
    
        gene_hits, gene_predictors, gene_pfixs = (numpy.array(x) for x in zip(*sorted(zip(gene_hits, gene_predictors, gene_pfixs), key=lambda pair: (pair[0]))))    
    
        early_predictors.sort()
        late_predictors.sort()
       
        rescaled_predictors = numpy.exp(numpy.fabs(numpy.log(predictors/mavg)))
        
        logit_mod = sm.Logit(responses, sm.add_constant(rescaled_predictors))
        logit_res = logit_mod.fit()

        sys.stdout.write("Logistic regression for %ss:\n" % metapopulation)
        sys.stdout.write("%s\n" % str(logit_res.summary()))
        sys.stdout.write("Params:\n")
        sys.stdout.write("%s\n" % str(logit_res.params))
        
        sys.stdout.write("Avg mutation multiplicity=%g, Avg fixed mutation multiplicity=%g\n" % (predictors.sum()/len(responses), (predictors*responses).sum()/responses.sum()))
        
        bins = numpy.logspace(0,2,10)
        bins[0] = 1e-02
        bin_values = numpy.exp(numpy.log(bins[1:])-(numpy.log(bins[2])-numpy.log(bins[1]))/2)
        
        
        bin_idxs = numpy.digitize(predictors,bins=bins)
 
        binned_num = numpy.array([0 for i in xrange(0,len(bins)-1)])
        binned_response = numpy.array([0 for i in xrange(0,len(bins)-1)])
    
        for i in xrange(0,len(predictors)):
            binned_num[bin_idxs[i]-1]+=1
            binned_response[bin_idxs[i]-1]+=responses[i]
        
        binned_pfix = binned_response*1.0/(binned_num+(binned_num==0))
    
        nonzero_bin_idxs = numpy.nonzero(binned_num)[0]
        lowest_nonzero_bin_idx = nonzero_bin_idxs[0]
        highest_nonzero_bin_idx = nonzero_bin_idxs[-1]
    
        theory_predictors = numpy.logspace(log10(bin_values[lowest_nonzero_bin_idx]/1.3),log10(bin_values[highest_nonzero_bin_idx]*1.3),100)
        theory_rescaled_predictors = numpy.exp(numpy.fabs(numpy.log(theory_predictors/mavg)))
        
        theory_response = 1/(1+numpy.exp(-(logit_res.params[0]+logit_res.params[1]*theory_rescaled_predictors)))
    
        # calculate percentile estimates for each bin using Bayesian posterior
        binned_pfix_upper = numpy.array([beta.ppf(0.84, binned_response[i]+syn_pfix, binned_num[i]-binned_response[i]+(1-syn_pfix), loc=0, scale=1) for i in xrange(0,len(binned_num))])
        binned_pfix_mid = numpy.array([beta.ppf(0.50, binned_response[i]+1, binned_num[i]-binned_response[i]+1, loc=0, scale=1) for i in xrange(0,len(binned_num))])
        binned_pfix_lower = numpy.array([beta.ppf(0.16, binned_response[i]+syn_pfix, binned_num[i]-binned_response[i]+(1-syn_pfix), loc=0, scale=1) for i in xrange(0,len(binned_num))])
    
        alpha=0.5
    
        for i in xrange(0,len(bin_values)):
            if binned_num[i]>0:
                pfix_axes[include_svs][metapopulation].bar(bins[i], binned_pfix[i]+0.05, (bins[i+1]-bins[i]),edgecolor='none',bottom=-0.05,linewidth=0, facecolor=colors[metapopulation],alpha=0.5)
                pfix_axes[include_svs][metapopulation].plot([bin_values[i],bin_values[i]], [binned_pfix_lower[i],binned_pfix_upper[i]],'-',color=colors[metapopulation],alpha=alpha,linewidth=0.5)
         
           
        pfix_axes[include_svs][metapopulation].plot(gene_predictors[-20:], gene_pfixs[-20:],'.',color=colors[metapopulation],alpha=alpha,markersize=2)
        
        sys.stdout.write("Dots for num_hits>=%g\n" % gene_hits[-20:].min())
        
        pfix_axes[include_svs][metapopulation].plot(theory_predictors, theory_response, '-', color=colors[metapopulation], alpha=alpha)
              
        line, = pfix_axes[include_svs][metapopulation].semilogx([1,100], [syn_pfix, syn_pfix],  '-', color=colors[metapopulation], linewidth=0.5)       
        
        line.set_dashes((3,2))
        pfix_axes[include_svs][metapopulation].set_ylim([0,1.05])
        pfix_axes[include_svs][metapopulation].set_xlim([1,100])
        pfix_axes[include_svs][metapopulation].set_yticklabels([])
        
        # "Done digitizing!"
    
        sys.stderr.write("Calculating null distribution...\n")
        null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)
        
        theory_ms = numpy.logspace(0,2,100)
        theory_survivals = null_survival_function(theory_ms)
        theory_survivals /= theory_survivals[0]
        
        sys.stderr.write("Done!\n")
    
        parallelism_axes[include_svs][metapopulation].step(predictors, (len(predictors)-numpy.arange(0,len(predictors)))*1.0/len(predictors), where='post',color=colors[metapopulation],alpha=alpha,label='All')
    
        parallelism_axes[include_svs][metapopulation].step(late_predictors, (len(late_predictors)-numpy.arange(0,len(late_predictors)))*1.0/len(late_predictors), where='post',color=colors[metapopulation],alpha=0.25,label='$>\mathrm{Median}(t)$')
    
        parallelism_axes[include_svs][metapopulation].step(early_predictors, (len(early_predictors)-numpy.arange(0,len(early_predictors)))*1.0/len(early_predictors), where='post',color=colors[metapopulation],alpha=0.75,label='$\leq \mathrm{Median}(t)$')
        
        parallelism_axes[include_svs][metapopulation].loglog(theory_ms,  theory_survivals,'-',alpha=alpha,linewidth=0.5,color='0.7')
    
        sys.stdout.write("Total mutations = %g\n" % len(predictors))
        sys.stdout.write("# m >= 2: %g\n" % (predictors>=2).sum())
        sys.stdout.write("%g of total\n" % ((predictors>=2).sum()*1.0/len(predictors)))
        sys.stdout.write("Expected # m >=2: %g\n" %  (theory_survivals[numpy.nonzero(theory_ms>=2)[0][0]]*len(predictors)))
        sys.stdout.write("Expected mavg=%g\n" % mavg)
        parallelism_axes[include_svs][metapopulation].set_ylim([1.0/len(predictors), 1.3])
        parallelism_axes[include_svs][metapopulation].set_xlim([1,100])
        
        if metapopulation=='mutator':
            legend_loc='center left'
        else:
            legend_loc='upper right'
            
        parallelism_axes[include_svs][metapopulation].legend(loc=legend_loc,frameon=False)

#####
#
# Tweak xlim,ylim and ticklabels
#
#####

###
#
# Fig 5
#
###

nonmutator_axis.set_xticks([i*10000 for i in xrange(0,7)]+[67500])
nonmutator_axis.set_xticklabels(["" for i in xrange(0,9)])
mutator_axis.set_xticks([i*10000 for i in xrange(0,7)]+[67500])
mutator_axis.set_xticklabels(['0']+['%d0k' % i for i in xrange(1,7)]+['Fixed'])
       
nonmutator_axis.legend(loc='center left',frameon=False)
vartype_pfix_axis.legend(loc='upper right',frameon=False)        

nonmutator_parallelism_axis.set_xlim([1,1e02])        
nonmutator_parallelism_axis.set_xticklabels([])
mutator_parallelism_axis.set_xlim([1,1e02])

nonmutator_pfix_axis.set_xlim([1,1e02])
nonmutator_pfix_axis.set_ylim([0,1.05])
nonmutator_pfix_axis.set_yticklabels([])

mutator_pfix_axis.set_xlim([1,1e02])
mutator_pfix_axis.set_ylim([0,1.05])
mutator_pfix_axis.set_yticklabels([])

nonmutator_axis.text(2000,1.9,figure_utils.get_panel_label('a'),fontsize=6,fontweight='bold')
mutator_axis.text(2000,32.5,figure_utils.get_panel_label('b'),fontsize=6,fontweight='bold')
vartype_pfix_axis.text(0,0.97,figure_utils.get_panel_label('c'),fontsize=6,fontweight='bold')
nonmutator_parallelism_axis.text(1.3,1.3,figure_utils.get_panel_label('d'),fontsize=6,fontweight='bold')
mutator_parallelism_axis.text(1.3,1.3,figure_utils.get_panel_label('e'),fontsize=6,fontweight='bold')
nonmutator_pfix_axis.text(1.5,0.97,figure_utils.get_panel_label('f'),fontsize=6,fontweight='bold')
mutator_pfix_axis.text(1.5,0.97,figure_utils.get_panel_label('g'),fontsize=6,fontweight='bold')

###
#
# No sv fig
#
###

nonmutator_parallelism_axis_nosv.set_xlim([1,1e02])        
nonmutator_parallelism_axis_nosv.set_xticklabels([])
mutator_parallelism_axis_nosv.set_xlim([1,1e02])
mutator_parallelism_axis_nosv.set_xticklabels([])

nonmutator_pfix_axis_nosv.set_xlim([1,1e02])
nonmutator_pfix_axis_nosv.set_ylim([0,1.05])

mutator_pfix_axis_nosv.set_xlim([1,1e02])
mutator_pfix_axis_nosv.set_ylim([0,1.05])
mutator_pfix_axis_nosv.set_yticklabels([])


# save figure to disk
fig5.savefig(parse_file.figure_directory+'fig5.pdf',bbox_inches='tight')  
nosv_fig.savefig(parse_file.figure_directory+'supplemental_parallelism_nosv.pdf',bbox_inches='tight')
#pylab.show()