from numpy.random import choice,multinomial
import numpy
from scipy.special import gammaln
from math import exp,log,fabs
from scipy.optimize import newton
import parse_file
import sys
from scipy.stats import poisson, binom
import stats_utils

def calculate_parallelism_statistics(convergence_matrix, allowed_populations=parse_file.complete_nonmutator_lines,Lmin=0):

    allowed_populations = set(allowed_populations)
    
    gene_statistics = {}
    
    # Now calculate gene counts
    Ltot = 0
    Ngenes = 0
    ntot = 0
    for gene_name in sorted(convergence_matrix.keys()):
        
        times = []
        
        L = max([convergence_matrix[gene_name]['length'],Lmin])
        n = 0
        
        num_pops = 0
        
        for population in allowed_populations:
            
            new_muts = len(convergence_matrix[gene_name]['mutations'][population])
            
            if new_muts > 0.5:
            
                num_pops += 1
            
                n += new_muts
            
                for t,l,lclade,f in convergence_matrix[gene_name]['mutations'][population]:
                    times.append(t)
            
        Ltot += L
        ntot += n
        
        gene_statistics[gene_name] = {}
        gene_statistics[gene_name]['length'] = L
        gene_statistics[gene_name]['observed'] = n
        if len(times) > 0:
            gene_statistics[gene_name]['median_time'] = numpy.median(times)
        else:
            gene_statistics[gene_name]['median_time'] = 0
            
        gene_statistics[gene_name]['nonzero_populations'] = num_pops
        
    Lavg = Ltot*1.0/len(convergence_matrix.keys())
    
    for gene_name in gene_statistics.keys():
    
        gene_statistics[gene_name]['expected'] = ntot*gene_statistics[gene_name]['length']/Ltot
        
        gene_statistics[gene_name]['multiplicity'] = gene_statistics[gene_name]['observed']*1.0/gene_statistics[gene_name]['length']*Lavg
        
        gene_statistics[gene_name]['g'] = log(gene_statistics[gene_name]['observed']*1.0/gene_statistics[gene_name]['expected']+(gene_statistics[gene_name]['observed']==0))
        
    return gene_statistics

def calculate_g_scores(gene_statistics, allowed_genes=None):
    # Calculate g score for each gene, i.e.
    # log(n_observed/n_expected)

    if allowed_genes==None:
        allowed_genes = set(gene_statistics.keys())
    else:
        allowed_genes = set(allowed_genes)
    
    Ltot = 0
    Lsig = 0
    nsig = 0
    ntot = 0
        
    for gene_name in gene_statistics.keys():
        
        Ltot += gene_statistics[gene_name]['length']
        ntot += gene_statistics[gene_name]['observed']
        
        if gene_name in allowed_genes:
            Lsig += gene_statistics[gene_name]['length']
            nsig += gene_statistics[gene_name]['observed']
            
    if fabs(nsig-ntot) < 0.5:
        normalization_factor = 1.0
    else:
        normalization_factor = (1-Lsig/Ltot)/(1-nsig*1.0/ntot)
        
    gene_g_scores = {}
    for gene_name in gene_statistics.keys():
    
        if gene_name in allowed_genes:
            
            if gene_statistics[gene_name]['observed'] < 0.5:
                gene_g_scores[gene_name] = 0
            else:
                gene_g_scores[gene_name] = log( gene_statistics[gene_name]['observed']*1.0/gene_statistics[gene_name]['expected']*normalization_factor)
        
        else:
            gene_g_scores[gene_name] = 0
            
    return gene_g_scores

def calculate_G_scores(gene_statistics, allowed_genes=None):
    # Calculates the G score for the whole gene, i.e.
    # n*g
    
    gene_g_scores = calculate_g_scores(gene_statistics,allowed_genes)
    
    gene_G_scores = {gene_name: gene_statistics[gene_name]['observed']*gene_g_scores[gene_name] for gene_name in gene_g_scores.keys()}
    
    return gene_G_scores        
        
def calculate_total_parallelism(gene_statistics, allowed_genes=None, num_bootstraps=10000):
    
    if allowed_genes==None:
        allowed_genes = gene_statistics.keys()
        
    Ls = []
    ns = []
    
    for gene_name in allowed_genes:
        
        Ls.append( gene_statistics[gene_name]['length'] )
        ns.append( gene_statistics[gene_name]['observed'] )

        
    Ls = numpy.array(Ls)
    ns = numpy.array(ns)
    
    Ltot = Ls.sum()
    ntot = ns.sum()
    ps = Ls*1.0/Ltot
    
    gs = ns*numpy.log(ns/(ntot*ps)+(ns==0))
    
    observed_G = gs.sum()/ns.sum()
    bootstrapped_Gs = []    
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_ns = multinomial(ntot,ps)
        bootstrapped_gs = bootstrapped_ns*numpy.log(bootstrapped_ns/(ntot*ps)+(bootstrapped_ns==0))
        bootstrapped_G = bootstrapped_gs.sum()/bootstrapped_ns.sum()
        
        bootstrapped_Gs.append(bootstrapped_G)
        
    bootstrapped_Gs = numpy.array(bootstrapped_Gs)
    
    pvalue = ((bootstrapped_Gs>=observed_G).sum()+1.0)/(len(bootstrapped_Gs)+1.0)
    return observed_G, pvalue

def calculate_parallelism_logpvalues(gene_statistics):

    gene_names = []
    Ls = []
    ns = []
    expected_ns = []
    
    for gene_name in gene_statistics.keys():        
        gene_names.append(gene_name)
        ns.append(gene_statistics[gene_name]['observed'])
        expected_ns.append(gene_statistics[gene_name]['expected'])
        
    ns = numpy.array(ns)
    expected_ns = numpy.array(expected_ns)

    logpvalues = stats_utils.calculate_poisson_log_survival(ns, expected_ns)
    
    return {gene_name: logp for gene_name, logp in zip(gene_names, logpvalues)}
        
def calculate_parallelism_qvalues(gene_statistics):

    gene_names = []
    Ls = []
    ns = []
    expected_ns = []
    
    for gene_name in gene_statistics.keys():        
        gene_names.append(gene_name)
        ns.append(gene_statistics[gene_name]['observed'])
        expected_ns.append(gene_statistics[gene_name]['expected'])
        
    ns = numpy.array(ns)
    expected_ns = numpy.array(expected_ns)
    ntot = ns.sum()
    ps = expected_ns/ntot
    ntots = ntot*numpy.ones_like(ps)

    pvalues = binom.sf(ns-0.5, ntots, ps)
    
    qvalues = stats_utils.calculate_qvalues(pvalues)
    
    qvalue_map = {gene_name: q for gene_name,q in zip(gene_names, qvalues)}
    pvalue_map = {gene_name: p for gene_name,p in zip(gene_names, pvalues)}
    
    return qvalue_map, pvalue_map 
    

class NullMultiplicitySurvivalFunction(object):
    # Null multiplicity distribution for mutations

    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls*1.0/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot)
    
    def __call__(self, m):
        
        lower_limits = numpy.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-2+0.1
        return (poisson.sf(lower_limits, self.expected_ns[None,:])*self.ps[None,:]).sum(axis=1)

class NullUniformGeneHitSurvivalFunction(object):
    # Null hit distribution for genes, 
    # but assume that all genes have same length
    
    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = numpy.ones_like(Ls)*1.0/len(Ls)
        self.expected_ns = self.ntot*self.ps
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot)
    
    def __call__(self, n):
        return (poisson.sf(n[:,None]-0.5, self.expected_ns[None,:])).sum(axis=1)


class NullGeneHitSurvivalFunction(object):
    # Null hit distribution for genes
    
    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls*1.0/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot)
    
    def __call__(self, n):
        return (poisson.sf(n[:,None]-0.5, self.expected_ns[None,:])).sum(axis=1)

        
class NullGeneMultiplicitySurvivalFunction(object):
    # Null multiplicity distribution for genes
    

    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls*1.0/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot)
    
    def __call__(self, m):
        
        lower_limits = numpy.ceil(m[:,None]*self.Ls[None,:]/self.Lavg)-1+0.1
        return (poisson.sf(lower_limits, self.expected_ns[None,:])).sum(axis=1)

class NullGeneGSurvivalFunction(object):

    # returns unnormalized survival function (raw number of genes)

    def __init__(self, Ls, ntot):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)*1.0
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot)
    
    def __call__(self, gs):
        
        # Do sum by hand
        ns = numpy.arange(0,400)*1.0
        gscores = ns[None,:]*numpy.log(ns[None,:]/self.expected_ns[:,None]+(ns[None,:]<0.5))
        
        logprobabilities = ns[None,:]*numpy.log(self.expected_ns)[:,None]-gammaln(ns+1)[None,:]-self.expected_ns[:,None]
        probabilities = numpy.exp(logprobabilities)
        survivals = numpy.array([ ((gscores>=g)*probabilities).sum() for g in gs])
        return survivals

class NullGeneLogpSurvivalFunction(object):
    # Null distribution of -log p for each gene
    
    def __init__(self, Ls, ntot,nmin=0):
        self.ntot = ntot
        self.Ls = numpy.array(Ls)*1.0
        self.Lavg = self.Ls.mean()
        self.ps = self.Ls/self.Ls.sum()
        self.expected_ns = self.ntot*self.ps
        self.nmin = nmin
    
    @classmethod
    def from_parallelism_statistics(cls, gene_parallelism_statistics,nmin=0):
        
        # calculate Ls
        Ls = []
        ntot = 0
        for gene_name in gene_parallelism_statistics.keys():
            Ls.append(gene_parallelism_statistics[gene_name]['length'])
            ntot += gene_parallelism_statistics[gene_name]['observed']
        
        return cls(Ls, ntot, nmin)
    
    def __call__(self, mlogps):
        
        # Do sum by hand
        ns = numpy.arange(0,400)*1.0
        
        logpvalues = stats_utils.calculate_poisson_log_survival(ns[None,:], self.expected_ns[:,None])
        
        logprobabilities = ns[None,:]*numpy.log(self.expected_ns)[:,None]-gammaln(ns+1)[None,:]-self.expected_ns[:,None]
        probabilities = numpy.exp(logprobabilities)
        survivals = numpy.array([ ((logpvalues>=mlogp)*(ns[None,:]>=self.nmin)*probabilities).sum() for mlogp in mlogps])
        return survivals


def calculate_scaled_missed_opportunities(pop_vector):
    num_hits = pop_vector.sum()*1.0
    num_nonzero_pops = (pop_vector>0.5).sum()*1.0
    num_pops = len(pop_vector)*1.0
    
    return (1.0-num_nonzero_pops/num_pops)*(num_hits-num_nonzero_pops)
    
def calculate_scaled_missed_opportunities_from_matrix(pop_matrix):
    ps = pop_matrix.sum(axis=0)*1.0
    ps /= ps.sum()
    
    p0s = (ps[None,:]*(pop_matrix<0.5)).sum(axis=1)
    
    num_hits = pop_matrix.sum(axis=1)*1.0
    num_nonzero_pops = (pop_matrix>0.5).sum(axis=1)*1.0
    num_pops = pop_matrix.shape[1]*1.0
    
    return p0s*(num_hits-num_nonzero_pops)
    #return (1.0-num_nonzero_pops/num_pops)*(num_hits-num_nonzero_pops)
    

def calculate_missed_opportunities(pop_vector):
    return min([len(pop_vector),pop_vector.sum()])-(pop_vector>0.5).sum()

def calculate_missed_opportunities_from_matrix(pop_matrix):
    return numpy.fmin(pop_matrix.sum(axis=1), pop_matrix.shape[1])-(pop_matrix>0.5).sum(axis=1)

def resample_population_matrix(population_matrix):

    ns = population_matrix.sum(axis=1)

    ps = population_matrix.sum(axis=0)*1.0
    ps /= ps.sum()
    
    bootstrapped_matrix = numpy.array([multinomial(n,ps) for n in ns])
    return bootstrapped_matrix
 
 
#######
# old stuff
    
def calculate_clade_grouping_loglikelihood(hit_table,ws):
    
    p = hit_table.shape[1]/2
    # loop over all groupings 
    
    total_associations = 2**(p-1)
    clade_associations = []
    for association_idx in xrange(0,total_associations):
        group_1 = [1,0]
        group_2 = [0,1]
        for bit in xrange(0,p-1):
            bit_set = not ((association_idx) & (1<<(bit)))
            if bit_set:
               group_1.extend([1,0])
               group_2.extend([0,1])
            else:
               group_1.extend([0,1])
               group_2.extend([1,0])
            
        clade_associations.append((numpy.array(group_1), numpy.array(group_2)))
    
    ps = ws*1.0/ws.sum()
    
    gene_hits = hit_table.sum(axis=1)
        
    loglikelihoods = []
        
    for group_1,group_2 in clade_associations:
    
        clade1_gene_hits = hit_table[:,group_1].sum(axis=1)
        clade2_gene_hits = hit_table[:,group_2].sum(axis=1)
        p1 = ps[group_1].sum()
        p2 = ps[group_2].sum()
    
        ell = (clade1_gene_hits*numpy.log((clade1_gene_hits+(clade1_gene_hits==0))*1.0/(gene_hits+(gene_hits==0))/p1)).sum() + (clade1_gene_hits*numpy.log((clade2_gene_hits+(clade2_gene_hits==0))*1.0/(gene_hits+(gene_hits==0))/p2)).sum() 
        loglikelihoods.append(ell)
        
    loglikelihoods = numpy.array(loglikelihoods)
    return loglikelihoods.max()

def clade_grouping_test(hit_table, num_bootstraps=1000):
    
    observed_loglikelihood = calculate_clade_grouping_loglikelihood(hit_table)
    
    ps = hit_table.sum(axis=0)*1.0/hit_table.sum()
    ntots = hit_table.sum(axis=1)
    
    bootstrapped_loglikelihoods = []
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_hit_table = []
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            bootstrapped_hit_table.append(ns)
            
        
        bootstrapped_hit_table = numpy.array(bootstrapped_hit_table)
        bootstrapped_loglikelihood = calculate_clade_grouping_loglikelihood(bootstrapped_hit_table)
        bootstrapped_loglikelihoods.append(bootstrapped_loglikelihood)
    
    bootstrapped_loglikelihoods = numpy.array(bootstrapped_loglikelihoods)
    
    pvalue = ((bootstrapped_loglikelihoods>=observed_loglikelihood).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)
    
    print observed_loglikelihood, bootstrapped_loglikelihoods.mean(), bootstrapped_loglikelihoods.std()*2, pvalue
    return observed_loglikelihood, pvalue
  

def calculate_log_hypergeometric_probability(ns,ws):
    
    ntot = ns.sum()
    wtot = ws.sum()
    
    ell = (gammaln(ws+1)-gammaln(ns+1)-gammaln(ws-ns+1)).sum()-gammaln(wtot+1)+gammaln(ntot+1)+gammaln(wtot-ntot+1)
    return (-ell)

def total_hypergeometric_test(hit_table,ws,num_bootstraps=1000):
    
    ps = ws*1.0/ws.sum()
    ntots = hit_table.sum(axis=1)
    
    observed_prob = 0
    for i in xrange(0,len(ntots)):
        observed_prob += calculate_log_hypergeometric_probability(hit_table[i],ws)
    
    bootstrapped_probs = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
        
        bootstrapped_prob = 0
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            bootstrapped_prob += calculate_log_hypergeometric_probability(ns,ws)
        
        bootstrapped_probs.append(bootstrapped_prob)
        
    bootstrapped_probs = numpy.array(bootstrapped_probs)
    
    pvalue = ((bootstrapped_probs>=observed_prob).sum()+1.0)/(len(bootstrapped_probs)+1.0)
    
    print observed_prob, bootstrapped_probs.mean(), bootstrapped_probs.std()*2, pvalue
    return observed_prob, pvalue
  
def total_zeros_multinomial_test(hit_table,ws,num_bootstraps=1000):
    
    ps = ws*1.0/ws.sum()
    ntots = hit_table.sum(axis=1)
    
    observed_prob = 0
    for i in xrange(0,len(ntots)):
        observed_prob += calculate_zeros_multinomial_distance(hit_table[i],ws)
    
    bootstrapped_probs = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
        
        bootstrapped_prob = 0
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            bootstrapped_prob += calculate_zeros_multinomial_distance(ns,ws)
        
        bootstrapped_probs.append(bootstrapped_prob)
        
    bootstrapped_probs = numpy.array(bootstrapped_probs)
    
    pvalue = ((bootstrapped_probs>=observed_prob).sum()+1.0)/(len(bootstrapped_probs)+1.0)
    
    return observed_prob, bootstrapped_probs.mean(), bootstrapped_probs.std()*2, pvalue


def calculate_multinomial_distance(ns,ws):
    ntot = ns.sum()*1.0
    wtot = ws.sum()*1.0
    distance = (ns*numpy.log((ns+(ns<0.5))*wtot/ntot/ws)).sum()
    return distance

def calculate_zeros_multinomial_distance(ns,ws):
    ntot = ns.sum()*1.0
    wtot = ws.sum()*1.0
    wtot_nonzero = ws[ns>0].sum()
    
    distance = ntot*log(wtot/wtot_nonzero)
    
    return distance/ntot

def calculate_entropy_distance(ns,ws):
    ntot = ns.sum()*1.0
    wtot = ws.sum()*1.0
    distance = ((ns/ntot)*numpy.log((ns+(ns<0.5))/ntot)).sum()-((ws/wtot)*numpy.log(ws/wtot)).sum()
    return distance

def zero_counts_test(hit_table,ws,num_bootstraps=1000):
    
    ps = ws*1.0/ws.sum()
    ntots = hit_table.sum(axis=1)
    
    observed_counts = [0 for i in xrange(0,len(ps))]
    
    for i in xrange(0,len(ntots)):
        observed_counts[(hit_table[i,:]==0).sum()]+=1
        
    observed_counts = numpy.array(observed_counts)*1.0
    
    bootstrapped_counts = [0 for i in xrange(0,len(ps))]
    
    for bootstrap_idx in xrange(0,num_bootstraps):    
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            bootstrapped_counts[(ns==0).sum()]+=1
            
    bootstrapped_counts = numpy.array(bootstrapped_counts)*1.0
    bootstrapped_counts = bootstrapped_counts/num_bootstraps
            
    return observed_counts, bootstrapped_counts
    
def excess_zeros_test(hit_table,ws,num_bootstraps=1000,alpha=0.05):
    
    ps = ws*1.0/ws.sum()
    ntots = hit_table.sum(axis=1)
    
    min_num_zeros = (len(ps)-ntots)*(ntots<len(ps))
    
    observed_excess_zeros = ((hit_table==0).sum(axis=1)-min_num_zeros).sum()
    
    bootstrapped_excess_zeros = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
        
        excess_zeros = 0.0
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            excess_zeros += ((ns==0).sum()-min_num_zeros[i])
        
        bootstrapped_excess_zeros.append(excess_zeros)
        
    bootstrapped_excess_zeros = numpy.array(sorted(bootstrapped_excess_zeros))
    
    percentile_idx = long(len(bootstrapped_excess_zeros)*alpha/2)
    
    null_mean = bootstrapped_excess_zeros.mean()
    null_std = bootstrapped_excess_zeros.std()
    null_lower = bootstrapped_excess_zeros[percentile_idx]
    null_upper = bootstrapped_excess_zeros[-percentile_idx]
    
    pvalue = ((bootstrapped_excess_zeros>=observed_excess_zeros).sum()+1.0)/(len(bootstrapped_excess_zeros)+1.0)
    
    #print observed_excess_zeros, bootstrapped_excess_zeros.mean(), bootstrapped_excess_zeros.std()*2, pvalue
    
    return observed_excess_zeros, null_mean, null_std, null_lower, null_upper, pvalue


def population_coincidence_test(hit_table,ws,num_bootstraps=1000):
    
    ps = ws*1.0/ws.sum()
    ntots = hit_table.sum(axis=1)
    
    #print "Population coincidence test"
    #print ps
    #print ntots
    
    coincidences = (hit_table>1.5).sum(axis=1)
    num_coincidences = (coincidences>0.5).sum()*1.0
    num_total = len(coincidences)*1.0
    
    observed_fraction = num_coincidences*1.0/num_total
    
    
    bootstrapped_fractions = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_num_coincidences = 0
        bootstrapped_num_total = 0
        for i in xrange(0,len(ntots)):
            ns = multinomial(ntots[i],ps)
            bootstrapped_num_total+=1
            if (ns>1.5).sum():
                bootstrapped_num_coincidences+=1
        bootstrapped_fractions.append(bootstrapped_num_coincidences*1.0/bootstrapped_num_total)
        
    bootstrapped_fractions = numpy.array(bootstrapped_fractions)
    pvalue = ((bootstrapped_fractions>=observed_fraction).sum()+1.0)/(len(bootstrapped_fractions)+1.0)
    print observed_fraction, bootstrapped_fractions.mean(), bootstrapped_fractions.std()*2, pvalue
    return observed_fraction, pvalue

def entropy_distance_test(ns,ws,num_bootstraps=1000):
    
    observed_distance = calculate_entropy_distance(ns,ws)
    
    bootstrapped_distances = []
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_ns = sample_null_hits(ns,ws)
        bootstrapped_distances.append( calculate_entropy_distance(bootstrapped_ns,ws) )
    
    bootstrapped_distances = numpy.array(bootstrapped_distances)
    pvalue = ((bootstrapped_distances>=observed_distance).sum()+1.0)/(len(bootstrapped_distances)+1.0)
    return observed_distance, pvalue
    
    
def multinomial_distance_test(ns,ws,num_bootstraps=1000):
    
    observed_distance = calculate_log_hypergeometric_probability(ns,ws)
    
    bootstrapped_distances = []
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_ns = sample_null_hits(ns,ws)
        bootstrapped_distances.append( calculate_log_hypergeometric_probability(bootstrapped_ns,ws) )
    
    bootstrapped_distances = numpy.array(bootstrapped_distances)
    pvalue = ((bootstrapped_distances>=observed_distance).sum()+1.0)/(len(bootstrapped_distances)+1.0)
    return observed_distance, pvalue
    
def sample_null_hits(ns, ws):
    
    N = ns.sum()
    W = ws.sum()*1.0
    
    # normalize ws
    normalized_ws = ws/W
    
    idxs = numpy.arange(0,len(ns))
    
    bootstrapped_ns = [0 for n in ns]
    
    chosen_idxs = choice(idxs,size=N,p=normalized_ws)
    
    for idx in chosen_idxs:
        bootstrapped_ns[idx] += 1
    
    bootstrapped_ns = numpy.array(bootstrapped_ns)
        
    return bootstrapped_ns

def zeros_multinomial_distance_test(ns,ws,num_bootstraps=1000):
    
    observed_distance = calculate_zeros_multinomial_distance(ns,ws)
    
    bootstrapped_distances = []
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_ns = sample_null_hits(ns,ws)
        bootstrapped_distances.append( calculate_zeros_multinomial_distance(bootstrapped_ns,ws) )
    
    bootstrapped_distances = numpy.array(bootstrapped_distances)
    pvalue = ((bootstrapped_distances>=observed_distance).sum()+1.0)/(len(bootstrapped_distances)+1.0)
    return observed_distance, bootstrapped_distances.mean(), bootstrapped_distances.std()*2, pvalue
    

def get_null_survival_function_from_total_genes(ns,ws):

    hs = numpy.arange(0,ns.max()+1)
    G = len(ns)*1.0
    mtot = (ns>0.5).sum()*1.0
    
    
    lam0 = ns.mean()/ws.mean()
    
    print lam0
    lam0 = newton(lambda x: numpy.exp(-x*ws).sum()-(G-mtot),lam0)
    print lam0
    
    phs = numpy.array([numpy.exp(h*numpy.log(lam0*ws)-gammaln(h+1)-lam0*ws).sum() for h in hs])
    
    survival_hs = numpy.array([phs[h:].sum() for h in hs]) 
    
    return hs, survival_hs

def get_null_scaled_survival_function_from_total_genes(ns,ws,wscale=1000.0):

    
    hs = numpy.arange(0,ns.max()+1)
    
    scaled_ns = ns*wscale/ws
    dh = scaled_ns[ns>0].min()
    theory_hs = dh*numpy.arange(0,numpy.ceil(scaled_ns.max()/dh)+1)
    
    print theory_hs[0:3]
    
    #theory_hs = numpy.linspace(0,scaled_ns.max(),100)

    #hs = numpy.arange(0,ns.max()+1)
    G = len(ns)*1.0
    mtot = (ns>0.5).sum()*1.0
    lam0 = ns.mean()/ws.mean()
    
    print lam0
    lam0 = newton(lambda x: numpy.exp(-x*ws).sum()-(G-mtot),lam0)
    print lam0
    
    phs = numpy.array([numpy.exp(h*numpy.log(lam0*ws)-gammaln(h+1)-lam0*ws).sum() for h in hs])
    
    individual_survivals = numpy.zeros((G,len(theory_hs)))
    for i in xrange(0,len(ws)):
        
        individual_phs = numpy.exp(hs*numpy.log(lam0*ws[i])-gammaln(hs+1)-lam0*ws[i])
        individual_effective_hs = hs*wscale/ws[i]
        
        individual_survivals[i,:] = numpy.array([individual_phs[individual_effective_hs>=h].sum() for h in theory_hs])
    
    
    survival_hs = individual_survivals.sum(axis=0)
    
    return theory_hs, survival_hs
     
def get_null_scaled_mutation_survival_function_from_total_genes(ns,ws,wscale=1000.0):

    
    hs = numpy.arange(0,ns.max()+1)
    
    scaled_ns = ns*wscale/ws
    dh = scaled_ns[ns>0].min()
    theory_hs = dh*numpy.arange(0,numpy.ceil(scaled_ns.max()/dh)+1)
    
    print theory_hs[0:3]
    
    #theory_hs = numpy.linspace(0,scaled_ns.max(),100)

    #hs = numpy.arange(0,ns.max()+1)
    G = len(ns)*1.0
    mtot = (ns>0.5).sum()*1.0
    lam0 = ns.mean()/ws.mean()
    
    print lam0
    lam0 = newton(lambda x: numpy.exp(-x*ws).sum()-(G-mtot),lam0)
    print lam0
    
    phs = numpy.array([numpy.exp(h*numpy.log(lam0*ws)-gammaln(h+1)-lam0*ws).sum() for h in hs])
    
    individual_survivals = numpy.zeros((G,len(theory_hs)))
    for i in xrange(0,len(ws)):
        
        individual_phs = numpy.exp(hs*numpy.log(lam0*ws[i])-gammaln(hs+1)-lam0*ws[i])
        individual_effective_hs = hs*wscale/ws[i]
        
        individual_survivals[i,:] = numpy.array([(hs*individual_phs)[individual_effective_hs>=h].sum() for h in theory_hs])
    
    
    survival_hs = individual_survivals.sum(axis=0)
    
    return theory_hs, survival_hs

def get_null_scaled_mutation_survival_function_per_site(ns, num_sites=4000000):
    lam0 = ns.sum()*1.0/num_sites
    hs = numpy.arange(0,ns.max()*2+1)*1.0
    phs = num_sites*numpy.exp(hs*numpy.log(lam0)-gammaln(hs+1)-lam0)

    survival_hs = numpy.array([(hs[hs>=h]*phs[hs>=h]).sum() for h in hs])
    return hs, survival_hs

def get_null_scaled_mutation_survival_function_from_total_mutations(ns,ws,normalization=None):

    if normalization==None:
        normalization = numpy.ones_like(ns)
    
    hs = numpy.arange(0,ns.max()+1)
    
    scaled_ns = ns/normalization
    dh = scaled_ns[ns>0].min()
    theory_hs = dh*numpy.arange(0,numpy.ceil(scaled_ns.max()/dh)+1)
    
    print theory_hs[0:3]
    
    #theory_hs = numpy.linspace(0,scaled_ns.max(),100)

    #hs = numpy.arange(0,ns.max()+1)
    G = len(ns)*1.0
    lam0 = ns.mean()/ws.mean()
    
    phs = numpy.array([numpy.exp(h*numpy.log(lam0*ws)-gammaln(h+1)-lam0*ws).sum() for h in hs])
    
    individual_survivals = numpy.zeros((G,len(theory_hs)))
    for i in xrange(0,len(ws)):
        
        individual_phs = numpy.exp(hs*numpy.log(lam0*ws[i])-gammaln(hs+1)-lam0*ws[i])
        individual_effective_hs = hs*normalization[i]
        
        individual_survivals[i,:] = numpy.array([(hs*individual_phs)[individual_effective_hs>=h].sum() for h in theory_hs])
    
    
    survival_hs = individual_survivals.sum(axis=0)
    
    return theory_hs, survival_hs
           

def get_null_survival_function_from_avg_hit(ns,ws):

    hs = numpy.arange(0,ns.max()+1)
    G = len(ns)
    
    observed_phs = numpy.array([((ns>(h-0.5))*(ns<(h+0.5))).sum() for h in hs])
    hstar = (hs*observed_phs).argmax()
    mstar = observed_phs[hstar]
    
    
    lam0 = ns.mean()/ws.mean()
    
    for iteration in xrange(0,100):
        lam0 = (mstar/numpy.exp(hstar*numpy.log(ws)-lam0*ws-gammaln(hstar+1)).sum())**(1.0/hstar)
    
         
    phs = numpy.array([numpy.exp(h*numpy.log(lam0*ws)-gammaln(h+1)-lam0*ws).sum() for h in hs])
    
    survival_hs = numpy.array([phs[h:].sum() for h in hs]) 
    
    return hs, survival_hs
    
def estimate_enriched_idxs(ns, ws, num_bootstraps=1000, FDR=0.5):

    # enriched relative to average

    idxs = numpy.arange(0,len(ns))
    
    N = ns.sum()*1.0
    W = ws.sum()*1.0
    
    nbars = N*ws/W 
    
    gs = ns*numpy.log((ns+(ns==0))/(nbars+(nbars==0)))
    
    sorted_gs = gs[ns>0]
    sorted_idxs = idxs[ns>0]
    
    sorted_gs, sorted_idxs = (numpy.array(x) for x in zip(*sorted(zip(sorted_gs, sorted_idxs), key=lambda pair: (pair[0]), reverse=True)))
    
    all_bootstrapped_gs = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        bootstrapped_ns = sample_null_hits(ns,ws)
        
        bootstrapped_gs = bootstrapped_ns*numpy.log((bootstrapped_ns+(bootstrapped_ns==0))/(nbars+(nbars==0)))
    
        sorted_bootstrapped_gs = bootstrapped_gs[bootstrapped_gs>0]
        all_bootstrapped_gs.extend(sorted_bootstrapped_gs)
        
    all_bootstrapped_gs = numpy.array(all_bootstrapped_gs)
    
    for i in xrange(0,len(sorted_idxs)):
        
        n_observed = (sorted_gs>=sorted_gs[i]).sum()
        n_bootstrapped = (all_bootstrapped_gs>=sorted_gs[i]).mean()*N
        
        
        if n_bootstrapped > (n_observed*FDR):
            return sorted_idxs[0:i], sorted_gs[i-1]
            
    return sorted_idxs, sorted_gs[-1]

def estimate_selected_idxs(ns, ws, num_bootstraps=1000, FDR=0.5):

    # selected relative to synonymous
    # idx = 0

    n_syn = ns[0]
    w_syn = ws[0]
    
    genic_ns = ns[1:]
    genic_ws = ws[1:]

    idxs = numpy.arange(1,len(ns))
    
    nbars = n_syn*genic_ws/w_syn 
    
    gs = genic_ns*numpy.log((genic_ns+(genic_ns==0))/(nbars+(nbars==0)))
    
    sorted_gs = gs[genic_ns>0]
    sorted_idxs = idxs[genic_ns>0]
    
    sorted_gs, sorted_idxs = (numpy.array(x) for x in zip(*sorted(zip(sorted_gs, sorted_idxs), key=lambda pair: (pair[0]), reverse=True)))
    
    all_bootstrapped_gs = []
    
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        bootstrapped_ns = sample_null_hits(ns,ws)
        
        bootstrapped_n_syn = bootstrapped_ns[0]
        bootstrapped_genic_ns = bootstrapped_ns[1:]
        bootstrapped_nbars = bootstrapped_n_syn * genic_ws / w_syn
        bootstrapped_gs = bootstrapped_genic_ns*numpy.log( (bootstrapped_genic_ns+(bootstrapped_genic_ns==0))/(bootstrapped_nbars+(bootstrapped_nbars==0)))
    
        sorted_bootstrapped_gs = bootstrapped_gs[bootstrapped_gs>0]
        all_bootstrapped_gs.extend(sorted_bootstrapped_gs)
        
    all_bootstrapped_gs = numpy.array(all_bootstrapped_gs)
    
    for i in xrange(0,len(sorted_idxs)):
        
        n_observed = (sorted_gs>=sorted_gs[i]).sum()
        n_bootstrapped = (all_bootstrapped_gs>=sorted_gs[i]).mean()*len(genic_ns)
        
        
        if n_bootstrapped > (n_observed*FDR):
            return sorted_idxs[0:i], sorted_gs[i-1]
            
    return sorted_idxs, sorted_gs[-1]
   
          
def estimate_mutation_spectrum_from_synonymous(ns, ws, selected_idxs=[]):

    rs = numpy.ones_like(ws)*1.0
    for idx in selected_idxs:
        rs[idx] = ns[idx]*ws[0]/ns[0]/ws[idx]
    
    ps = rs*ws
    ps /= ps[1:].sum()
    gs = ns*numpy.log(rs)
    
    return ps,rs,gs

def estimate_deltas_from_synonymous(nws, nms, ws, selected_idxs=[]):
    # deltas are first
    deltaws = numpy.zeros_like(ws)*1.0
    deltam1s = numpy.zeros_like(ws)*1.0
    deltam2s = numpy.zeros_like(ws)*1.0
    
    for i in selected_idxs:
        deltaws[i] = (nws[i]*1.0/ws[i]-nws[0]*1.0/ws[0])/(nws[0]*1.0/ws[0])
        deltam1s[i] = (nms[i]*1.0/ws[i]-nms[0]*1.0/ws[0])/(nms[0]*1.0/ws[0])
        deltam2s[i] = (nms[i]*1.0/ws[i]-nms[0]*1.0/ws[0])/(nws[0]*1.0/ws[0])
    
    return deltaws, deltam1s, deltam2s
    
def estimate_mutation_spectrum_from_relative(ns,ws,selected_idxs=[]):
    # estimate the selected ones
    
    N = ns.sum()*1.0
    W = ws.sum()*1.0
    
    nrest = N
    wrest = w
    for idx in selected_idxs:
        nrest -= ns[idx]
        wrest -= ws[idx]
    
    rs = numpy.ones_like(ws)*1.0
    for idx in selected_idxs:
        rs[idx] = ns[idx]*wrest/nrest/ws[idx]
    
    ps = rs*ws
    ps /= ps.sum()
    gs = ns*numpy.log(rs)
    
    return ps,rs,gs


def calculate_LRT_statistic(n1s,n2s):
    
    ntots = n1s+n2s
    
    ntot = ntots.sum()*1.0
    n1tot = n1s.sum()*1.0
    n2tot = n2s.sum()*1.0
    
    if n1tot < 0.5:
        return 0
        
    if n2tot < 0.5:
        return 0
    
    pbars = ntots/ntot
    p1s = n1s/n1tot
    p2s = n2s/n2tot
    
    #print ntot,n1tot,n2tot
    
    LRT = 0
    # counting factors:
    #LRT += gammaln(n1tot+1)+gammaln(n2tot+1)-gammaln(ntot+1)
    #LRT += (gammaln(ntots*(ntots>0.5)+1)-gammaln(n1s*(n1s>0.5)+1)-gammaln(n2s*(n2s>0.5)+1)).sum()
    # probability factors
    LRT += -1*(ntots*(ntots>0.5)*numpy.log(pbars+(ntots<0.5))).sum()
    LRT += (n1s*(n1s>0.5)*numpy.log(p1s+(n1s<0.5))).sum()
    LRT += (n2s*(n2s>0.5)*numpy.log(p2s+(n2s<0.5))).sum()
    #print LRT
    return LRT 
    
def estimate_mutation_spectrum_from_synonymous_EM(all_ns, all_ws, deltas, r=1, num_iterations = 100):

    K = len(deltas)
    G = len(all_ns)-1
    
    w_syn = all_ws[0]*1.0
    n_syn = all_ns[0]*1.0
    
    ns = all_ns[1:]*1.0
    ws = all_ws[1:]*1.0
    
    ntot = ns.sum()*1.0
    wtot = ws.sum()*1.0
    
    # initial estimates for lam0 and pks
    
    lam0 = n_syn/w_syn
    rks = r+deltas
    
    pks = numpy.ones_like(rks)*1.0
    pks[0] = G
    pks/pks.sum()
    logrks = numpy.log(rks) 
    
    wbyrs = numpy.outer(ws,rks)
      
    for iteration in xrange(0,num_iterations): 
        
        # calculate qik (in logspace first)
        
        logqiks = numpy.zeros((G,K))
        logqiks = numpy.outer(ns, logrks)-lam0*wbyrs+numpy.outer(numpy.ones(G),numpy.log(pks))
    
        logqikmaxs = numpy.outer(logqiks.max(axis=1),numpy.ones(K))
        # subtract max
        qiks = numpy.exp(logqiks-logqikmaxs)
        # normalize
        qiks /= numpy.outer(qiks.sum(axis=1),numpy.ones(K))
        
        # calculate new estimate of lam0 and pk
        
        lam0 = (ntot+n_syn)/((qiks*wbyrs).sum()+w_syn)
        
        pks = qiks.sum(axis=0)/G
        
        if pks.argmax()!=zero_idx:
            pks = numpy.roll(pks,zero_idx-pks.argmax())
        
        # repeat
    
    return qiks, lam0, pks 
    
def estimate_mutation_spectrum_EM(ns, ws, deltas, r=1, num_iterations = 100):

    K = len(deltas)
    G = len(ns)
    
    zero_idx = numpy.fabs(deltas).argmin()
    
    ntot = ns.sum()*1.0
    wtot = ws.sum()*1.0
    
    # initial estimates for lam0 and pks
    
    lam0 = ntot/wtot
    rks = r+deltas
    
    pks = numpy.ones_like(rks)*1.0
    pks[zero_idx] = G
    pks/pks.sum()
    logrks = numpy.log(rks) 
    
    wbyrs = numpy.outer(ws,rks)
      
    for iteration in xrange(0,num_iterations): 
        
        # calculate qik (in logspace first)
        
        logqiks = numpy.zeros((G,K))
        logqiks = numpy.outer(ns, logrks)-lam0*wbyrs+numpy.outer(numpy.ones(G),numpy.log(pks))
    
        logqikmaxs = numpy.outer(logqiks.max(axis=1),numpy.ones(K))
        # subtract max
        qiks = numpy.exp(logqiks-logqikmaxs)
        # normalize
        qiks /= numpy.outer(qiks.sum(axis=1),numpy.ones(K))
        
        # calculate new estimate of lam0 and pk
        
        lam0 = (ntot)/((qiks*wbyrs).sum())
        
        pks = qiks.sum(axis=0)/G
        
        if pks.argmax()!=zero_idx:
            pks = numpy.roll(pks,zero_idx-pks.argmax())
        
        # repeat
    
    return qiks, lam0, pks 