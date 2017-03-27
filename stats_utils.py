import numpy
from numpy.random import uniform
from scipy.special import digamma,gammaln
from math import log
from scipy.stats import ks_2samp, kstest, poisson

tail_probability=0.25

def calculate_empirical_pvalue(observed_x, bootstrapped_xs):
    pvalue = ((bootstrapped_xs>=observed_x).sum()+1.0)/(len(bootstrapped_xs)+1.0)
    return pvalue

def calculate_qvalues(pvalues):

    # Uses the formula
    #
    # Qi = min_{Q>Pi} { Q*(sum_j 1)/(sum_j theta(Q-Pj)) } 
    #
    # Qi = max_{P_(k)>Pi} 

    # calculate q-values
    qvalues = []
    
    sorted_pvalues = numpy.array(sorted(pvalues))

    Ntot = len(sorted_pvalues)
    Nless = numpy.array([(sorted_pvalues<=p).sum() for p in sorted_pvalues])

    for p in pvalues:

        min_q = 1e06
            
        for j in reversed(xrange(0,Ntot)):
             
            if sorted_pvalues[j]<p:
                break
            else:
                new_q = Ntot*sorted_pvalues[j]*1.0/Nless[j]
                if new_q<min_q:
                    min_q = new_q
                    
        qvalues.append(min_q)
    
    qvalues = numpy.array(qvalues)
    
    return qvalues

def calculate_poisson_log_survival(ns, expected_ns):

    survivals = poisson.sf(ns-0.1, expected_ns)
    
    logsurvivals = numpy.zeros_like(survivals)
    logsurvivals[survivals>1e-20] = -numpy.log(survivals[survivals>1e-20])
    logsurvivals[survivals<=1e-20] = (-ns*numpy.log(ns/expected_ns+(ns==0))+ns-expected_ns)[survivals<=1e-20]

    return logsurvivals
    

def calculate_median(ts):
    # requires sorted ts
    n = len(ts)
    if (n%2)==0:
        median = (ts[n/2-1]+ts[n/2])/2.0
    else:
        median = ts[(n-1)/2]*1.0
        
    return median

def calculate_theoretical_ks_distance(xs, cdf):
    # 1-sample KS distance
    
    n = len(xs)*1.0
    return kstest(xs,cdf)[0]*(n)**0.5

def calculate_ks_distance(t1s,t2s):
    # 2-sample KS distance

    m = len(t1s)*1.0
    n = len(t2s)*1.0
    
    scipy_ks = ks_2samp(t1s, t2s)[0]*(n*m/(n+m))**0.5
    
    return scipy_ks
    
    # Old version, by hand (required sorted ts)
    # combine t1s and t2s 
    #ts = set(t1s)
    #ts.update(t2s)
    #n = len(t1s)
    #m = len(t2s)
    
    #D = numpy.fabs([(t1s<=t).mean()-(t2s<=t).mean() for t in ts]).max()
    #own_ks =  D*(m*n*1.0/(n+m))**0.5
    
    #print D, own_ks, scipy_ks
    
    #return scipy_ks

####
#
# Calculates unnormalized survival functions (i.e. # of observations >= x)
# from numpy vector of observations
#
####
def calculate_unnormalized_survival_from_vector(xs, min_x=None, max_x=None, min_p=1e-10):

    if min_x==None:
        min_x = xs.min()-1
    
    if max_x==None:
        max_x = xs.max()+1
        
    unique_xs = set(xs)
    unique_xs.add(min_x)
    unique_xs.add(max_x)
    
    xvalues = []
    num_observations = []
    
    for x in sorted(unique_xs):
        xvalues.append(x)
        num_observations.append( (xs>=x).sum() )
    
    # So that we can plot CDF, SF on log scale
    num_observations[0] -= min_p
    num_observations[1] -= min_p
    num_observations[-1] += min_p    
    
    return numpy.array(xvalues), numpy.array(num_observations)
 

####
#
# Calculates percentiles from sorted vector
#
####
def calculate_percentile_from_sorted_vector(sorted_xs, p):  
    return sorted_xs[long(len(sorted_xs)*p)]



def calculate_anderson_darling_distance(t1s,t2s):
    # combine t1s and t2s 
    ts = set(t1s)
    ts.update(t2s)
    ts = sorted(ts)
    
    n = len(t1s)*1.0
    m = len(t2s)*1.0
    
    Ns = numpy.array([(t1s<=t).sum() for t in ts])[:-1]
    Ms = numpy.array([(t2s<=t).sum() for t in ts])[:-1]
    
    Hs = (Ns+Ms)/(n+m) 
    dHs = numpy.array([(t1s==t).sum()+(t2s==t).sum() for t in ts])[:-1]
    
    Ds = Ns/n-Ms/m
    
    A = (numpy.square(Ds)/(Hs*(1-Hs))).sum()*n*m/(n+m)
    return A
    
# requires sorted ts
def calculate_variation_statistic(ts):
    
    n = len(ts)
    
    if n==1:
        return 0
    
    median = calculate_median(ts)
              
    tail_idx = long(n*tail_probability)
    range = ts[-1-tail_idx]-ts[tail_idx]+250   
  
    variation_statistic = range*1.0/median

    if n==2:
        variation_statistic /= min([2,2*(60000/median-1)])

    return variation_statistic 

def calculate_entropy(ts):
    
    n = len(ts)
    delta = 250.0
    
    perturbed_ts = ts+uniform(low=-250,high=250,size=n)
    perturbed_ts.sort()
    
    
    if n==1:
        return 0
    else:
        return (numpy.log(numpy.diff(perturbed_ts))).sum()/(n-1)+digamma(1)-digamma(n)


def estimate_mutation_spectrum(num_hits, weights):
    
    W = weights.sum()*1.0
    N = num_hits.sum()*1.0
    
    expected_num_hits = weights*N/W
    
    gs = num_hits*numpy.log(num_hits/expected_num_hits)
    
    sorted_num_hits = numpy.array(sorted(num_hits,reverse=True))
    G = num_genes*1.0
    ntot = sorted_num_hits.sum()
    hs = set(sorted_num_hits)
    hs = sorted(hs,reverse=True)
    # number of hits
    loglikelihoods = []
    fks = []
    fss = []
    
    for h in hs[:-1]:
        idxs = sorted_num_hits>=h
        g = (idxs).sum()
        nk = sorted_num_hits[idxs]
        fk = nk*1.0/nk.sum()
        ns = nk.sum()*1.0
        fs = 1.0+(G*ns/g/ntot-1)/(1-ns/ntot)
        po = 1.0/(G+g*(fs-1))
        
        loglikelihood = ns*log(g*fs)+ntot*log(po)+(nk*numpy.log(fk+(nk==0))).sum()+gammaln(g+1)+gammaln(G-g+1)+gammaln(g)-g*log(2)
        
        loglikelihoods.append(loglikelihood)
        fks.append(fk)
        fss.append(fs)
    
    loglikelihoods = numpy.array(loglikelihoods)
    h_idx = loglikelihoods.argmax()    
    print hs[h_idx],fss[h_idx]
    print loglikelihoods
        
            
        