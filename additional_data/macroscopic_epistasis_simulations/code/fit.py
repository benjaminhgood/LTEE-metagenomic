import numpy
import pylab
from numpy.random import multinomial, normal
from math import log,exp,sqrt
from scipy.optimize import fmin

def linear_trajectory(ts,m,b):
    return m*ts+b
    
def fit_linear_trajectory(ts,xs, mrange=None):
    m = ((xs*ts).mean()-xs.mean()*ts.mean())/((ts*ts).mean()-ts.mean()*ts.mean())
    
    if mrange != None:
        mmin = mrange[0]
        mmax = mrange[1]
        
        if m < mmin:
            m = mmin
        if m > mmax: 
            m = mmax
            
    b = (xs-m*ts).mean()
    return m,b
    
def powerlaw_trajectory(ts,loga,logb):
    result = numpy.log(1+exp(logb)*ts)*exp(loga)
    return result
    
def frank_mut_trajectory(ts,logXc,logv0,logUn):
    result = numpy.log(1+exp(logv0-logXc)*ts)+ts*exp(logUn)
    return result

def wiser_mut_trajectory(ts,logXc,logv0,logsf,logUn):
    result = 2*exp(logXc-logsf)*(numpy.sqrt(1+exp(logv0-logXc)*ts)-1)+ts*exp(logUn)
    return result
    
def finite_mut_trajectory(ts, logR0):
    result = ts*exp(logR0)
    return result
    
def coupon_trajectory(ts,loga,logb):
    result = exp(loga)*(1.0-1.0/numpy.square(1+exp(logb)*ts))
    return result
  
def sergey_trajectory(ts,loga,logb):
    result = exp(loga+logb)*ts/(exp(loga)+exp(logb)*ts)  
    return result
    
def fit_powerlaw_trajectory(ts,xs):
    loga0 = 0
    logb0 = log((exp(xs[-1])-1)/ts[-1])
    xmin = fmin(lambda x: numpy.square(xs-powerlaw_trajectory(ts,x[0],x[1])).sum(),numpy.array([loga0,logb0]))
    a = exp(xmin[0])
    b = exp(xmin[1])
    return a,b

def fit_sergey_trajectory(ts,xs):
    loga0 = log(5*xs[-1])
    logb0 = log(xs[-1]*exp(loga0)/(exp(loga0)-xs[-1])/ts[-1])
    xmin = fmin(lambda x: numpy.square(xs-sergey_trajectory(ts,x[0],x[1])).sum(),numpy.array([loga0,logb0]))
    a = exp(xmin[0])
    b = exp(xmin[1])
    return a,b

def fit_coupon_trajectory(ts,xs):
    loga0 = log(5*xs[-1])
    logb0 = log(xs[-1]*exp(loga0)/(exp(loga0)-xs[-1])/ts[-1])
    xmin = fmin(lambda x: numpy.square(xs-coupon_trajectory(ts,x[0],x[1])).sum(),numpy.array([loga0,logb0]))
    a = exp(xmin[0])
    b = exp(xmin[1])
    return a,b

def fit_wiser_mut_trajectory(ts,ms,logXc,logv0,Unmax=None):
    fs = 2*exp(logXc)*(numpy.sqrt(1+exp(logv0-logXc)*ts)-1)
    sf = (numpy.dot(fs,fs)*numpy.dot(ts,ts)-(numpy.dot(fs,ts))**2)/(numpy.dot(ms,fs)*numpy.dot(ts,ts)-numpy.dot(ms,ts)*numpy.dot(fs,ts))
    Un = (numpy.dot(ms,ts)-numpy.dot(fs,ts)/sf)/numpy.dot(ts,ts)
    
    if Unmax != None and Un > Unmax:
        Un = Unmax
        sf = numpy.dot(fs,fs)/(numpy.dot(ms,fs)-Un*numpy.dot(fs,ts))
    return sf,Un
    
def fit_finite_mut_trajectory(ts,ms,logXc,logv0):
    R0 = numpy.dot(ms, ts)/numpy.dot(ts,ts)
    return R0

def fit_neutral_mutation_rate(times, theory_muts, observed_muts):
    return numpy.dot(observed_muts-theory_muts, times)/numpy.dot(times,times)

# simulated trajectories is matrix (rows = independent runs, cols = timepoints)


#### ignore ###########################

def calculate_simulation_mse_pvalue(simulated_trajectories, observed_trajectory, sigma_measurement,num_bootstraps = 1000):
    avg_simulated_trajectory = simulated_trajectories.mean(axis=0)
    observed_mse = numpy.square(observed_trajectory-avg_simulated_trajectory).sum()
    bootstrapped_mses = []
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(6.0, simulation_weights)/6.0
        bootstrapped_observed_trajectory = numpy.dot(simulation_idxs, simulated_trajectories)+normal(0,sigma_measurement,size=len(observed_trajectory))
        bootstrapped_mses.append( numpy.square(bootstrapped_observed_trajectory-avg_simulated_trajectory).sum())
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    return observed_mse, (bootstrapped_mses > observed_mse).sum()*1.0/len(bootstrapped_mses)

def calculate_simulation_linear_pvalue(times, simulated_trajectories, observed_trajectory, sigma_measurement, num_bootstraps = 1000):
    m,b = fit_linear_trajectory(times,observed_trajectory)
    linear_fit = linear_trajectory(times,m,b)
    observed_mse = numpy.square(observed_trajectory-linear_fit).sum()
    bootstrapped_mses = []
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(6.0, simulation_weights)/6.0
        bootstrapped_observed_trajectory = numpy.dot(simulation_idxs, simulated_trajectories)+normal(0,sigma_measurement,size=len(observed_trajectory))
        m,b = fit_linear_trajectory(times,bootstrapped_observed_trajectory)
        linear_fit = linear_trajectory(times,m,b)
        bootstrapped_mses.append( numpy.square(bootstrapped_observed_trajectory-linear_fit).sum())
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    return observed_mse, (bootstrapped_mses > observed_mse).sum()*1.0/len(bootstrapped_mses)

def calculate_simulation_likelihood(simulated_trajectories, observed_trajectory, sigma_measurement, num_bootstraps = 1000):
    bootstrapped_mses = []
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(6.0, simulation_weights)/6.0
        bootstrapped_observed_trajectory = numpy.dot(simulation_idxs, simulated_trajectories)
        bootstrapped_mses.append( numpy.square(bootstrapped_observed_trajectory-observed_trajectory).sum())
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    individual_log_likelihoods = -bootstrapped_mses/2/sigma_measurement/sigma_measurement
    max_log_likelihood = individual_log_likelihoods.max()
    individual_log_likelihoods -= max_log_likelihood
    log_likelihood = max_log_likelihood + log( (numpy.exp(individual_log_likelihoods)).mean() ) 
    #print max_log_likelihood, log_likelihood
    return log_likelihood

################ end ignore ######################

def calculate_simulation_loglikelihood_pvalue(simulated_trajectories, observed_trajectory, sigma_measurement, simulated_dxs, observed_dx, sigma_dx_measurement, intercept=None, num_bootstraps = 10000):

    avg_simulated_trajectory = simulated_trajectories.mean(axis=0)
    avg_simulated_dx = simulated_dxs.mean()
    rescaled_sigma = sigma_measurement/sqrt(1.0*len(observed_trajectory))
    
    observed_mse = numpy.square(observed_trajectory-avg_simulated_trajectory).mean()
    observed_me = (observed_trajectory-avg_simulated_trajectory).mean()
    observed_dx_mse = (observed_dx-avg_simulated_dx)**2
    
    individual_mses = []
    individual_mes = []
    individual_dx_mses = []
    
    
    bootstrapped_mses = []
    bootstrapped_dx_mses = []
    
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(6.0, simulation_weights)/6.0
        trajectory = numpy.dot(simulation_idxs, simulated_trajectories)
        dx = numpy.dot(simulation_idxs, simulated_dxs)
        
        individual_mses.append( numpy.square(observed_trajectory-trajectory).mean() )
        individual_mes.append( (observed_trajectory-trajectory).mean() )
        individual_dx_mses.append( (observed_dx-dx)**2 )
        
        bootstrapped_mses.append( numpy.square(trajectory+normal(0,sigma_measurement,size=len(observed_trajectory))-avg_simulated_trajectory).mean() )
        
        bootstrapped_dx_mses.append( (dx + normal(0,sigma_dx_measurement))**2 )
    
    individual_mses = numpy.array(individual_mses)
    individual_mes = numpy.array(individual_mes)
    individual_dx_mses = numpy.array(individual_dx_mses)
    
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    bootstrapped_dx_mses = numpy.array(bootstrapped_dx_mses)
        
    if intercept == None:
        b0 = observed_me
        bs = numpy.linspace(0,0.3,1000)
    else:
        bs = numpy.array([intercept])
        
    def loglikelihood(b):
        individual_loglikelihoods = -(individual_mses-2*b*individual_mes+b*b)/2/rescaled_sigma/rescaled_sigma - individual_dx_mses/2/sigma_dx_measurement/sigma_dx_measurement
             
        return individual_loglikelihoods.max()+log(numpy.exp(individual_loglikelihoods-individual_loglikelihoods.max()).mean())
        
    loglikelihoods = numpy.array([loglikelihood(b) for b in bs])
    
    b = bs[loglikelihoods.argmax()]
    loglikelihood = loglikelihoods.max()
    observed_chisquared = (observed_mse - 2 * b * observed_me + b*b)/2/rescaled_sigma/rescaled_sigma + observed_dx_mse/2/sigma_dx_measurement/sigma_dx_measurement
    
    bootstrapped_chisquareds = bootstrapped_mses/2/rescaled_sigma/rescaled_sigma+bootstrapped_dx_mses/2/sigma_dx_measurement/sigma_dx_measurement
    
    pvalue = (bootstrapped_chisquareds > observed_chisquared).sum()*1.0/len(bootstrapped_chisquareds)
    
    return loglikelihood, pvalue, b   

def calculate_simulation_linear_pvalue(times, simulated_trajectories, observed_trajectory, sigma_measurement, num_bootstraps = 1000):
    m,b = fit_linear_trajectory(times,observed_trajectory)
    linear_fit = linear_trajectory(times,m,b)
    observed_mse = numpy.square(observed_trajectory-linear_fit).sum()
    bootstrapped_mses = []
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(6.0, simulation_weights)/6.0
        bootstrapped_observed_trajectory = numpy.dot(simulation_idxs, simulated_trajectories)+normal(0,sigma_measurement,size=len(observed_trajectory))
        m,b = fit_linear_trajectory(times,bootstrapped_observed_trajectory)
        linear_fit = linear_trajectory(times,m,b)
        bootstrapped_mses.append( numpy.square(bootstrapped_observed_trajectory-linear_fit).sum())
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    return observed_mse, (bootstrapped_mses > observed_mse).sum()*1.0/len(bootstrapped_mses)
    
def calculate_simulation_mutation_pvalue(simulated_trajectories, observed_trajectories, num_bootstraps=1000):
    avg_simulated_trajectory = simulated_trajectories.mean(axis=0)
    observed_mse = numpy.square(observed_trajectory/avg_simulated_trajectory-1).sum()
    bootstrapped_mses = []
    simulation_weights = numpy.ones(simulated_trajectories.shape[0])*1.0
    simulation_weights /= simulation_weights.sum()
    for i in xrange(0,num_bootstraps):
        simulation_idxs = multinomial(1.0, simulation_weights)/1.0
        bootstrapped_observed_trajectory = numpy.dot(simulation_idxs, simulated_trajectories)+normal(0,sigma_measurement,size=len(observed_trajectory))
        bootstrapped_mses.append( numpy.square(bootstrapped_observed_trajectory-avg_simulated_trajectory).sum())
    bootstrapped_mses = numpy.array(bootstrapped_mses)
    return observed_mse, (bootstrapped_mses > observed_mse).sum()*1.0/len(bootstrapped_mses)
    
def empirical_cdf(values,vmin=None,vmax=None):
    values = numpy.array(values)
    L = len(values)
    sorted_values = numpy.sort(values)
    values = numpy.zeros(2*L)
    values[0::2] = sorted_values
    values[1::2] = sorted_values
    
    cdf = numpy.zeros(2*L) 
    print len(cdf[1::2]), len(numpy.arange(1,L+1))
    cdf[0::2] = numpy.arange(1,L+1)*1.0-1
    cdf[1::2] = numpy.arange(1,L+1)*1.0
    cdf /= cdf[-1]
    if vmin!=None:
        values = numpy.hstack([numpy.array([vmin]), values])
        cdf = numpy.hstack([numpy.array([0]), cdf])
    if vmax!=None:
        values = numpy.hstack([values, numpy.array([vmax])])
        cdf = numpy.hstack([cdf, numpy.array([1])])
    return values, cdf
