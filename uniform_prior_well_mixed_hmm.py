import timecourse_utils
import parse_file
from math import log,exp
import numpy
from scipy.special import gammaln as loggamma
from scipy.special import betainc as incomplete_beta
import sys
from scipy.optimize import fmin   
from scipy.interpolate import interp1d
from scipy.stats import nbinom     

states = parse_file.well_mixed_hmm_states 

def infer_hmm(A, D, num_iterations=10):
    
    # run baum welch algorithm
     
    # fixed params
    epsF = 0.025 # max error rate from 100% for fixed mutations
    epsE = 0.025 # max error rate from 0% for extinct or unborn mutations
    min_error = 1e-03 # error rate can't be any lower than this (~ Illumina error rate)
    
    # transition rates for markov chain 
    birth_rate = 1e-02
    recurrence_rate = 1e-06
    extinction_rate = 0.5
    extinction_probability = 0.5
    
    M = A.shape[0] # mutations
    T = A.shape[1]     # timepoints
    ts = numpy.arange(0,T)
    L = len(states.keys())
    

    R = D-A
    safe_D = D+(D<0.5) # remove zero depths so that logs and stuff are safe. no influence on results
    safe_f = A/safe_D # remove zero freqs so that logs and stuff are safe. no influence on results
    safe_f = numpy.clip(safe_f,1e-03,1-1e-03)
    
    p0 = numpy.ones_like(safe_D)*1e-02
    p1 = numpy.ones_like(safe_D)*1e-02
    
    Q = numpy.zeros((T,L,L)) # transition matrix
    F = numpy.zeros((M,T,L))*1.0 # forward probability table
    B = numpy.zeros((M,T,L))*1.0 # backward probability table
    
    max_likelihood_table = numpy.zeros((M,T,L))*1.0
    state_pointer_table = numpy.zeros((M,T,L))
    
    Pstate = numpy.zeros((M,T,L))*1.0 # decoding probabilities for each state
        
    for current_iteration in xrange(1,num_iterations+1):
    
        sys.stderr.write("Iteration %d\n" % current_iteration)
        
        # calculate Qs and Vs
        #sys.stderr.write("Calculating Q,V...\n")
        
        Q = numpy.zeros((T, L, L))
        
        for t in xrange(0,T):
            
            # Unborn
            Q[t,states['A'],:] = numpy.array([1-birth_rate, recurrence_rate, 0, 0])
            
            # Extinct
            Q[t,states['E'],:] = numpy.array([0, 1-recurrence_rate, 0, extinction_rate*extinction_probability])
            
            # Fixed
            Q[t,states['F'],:] = numpy.array([0,0,1,extinction_rate*(1-extinction_probability)])
            
            # Polymorphic
            Q[t,states['P'],:] = numpy.array([birth_rate, 0, 0, 1-extinction_rate])
            
        
        logV = numpy.zeros((M,T,L))*1.0
        V = numpy.zeros_like(logV)
        
        #sys.stderr.write("Calculating emission probabilities...\n")
        
        # Unborn
        safe_unborn_beta = numpy.clip(incomplete_beta(A+1,R+1,2*p0),1e-300,1)
        logV[:,:,states['A']] = (D>0.5)*( numpy.log(safe_unborn_beta) - numpy.log(2*p0) ) 
        
        # Extinct 
        safe_extinct_beta = numpy.clip(incomplete_beta(A+1,R+1,2*p0),1e-300,1)
        logV[:,:,states['E']] = (D>0.5)*( numpy.log(safe_extinct_beta) - numpy.log(2*p0) ) 
        
        
        # Fixed 
        safe_fixed_beta = numpy.clip(incomplete_beta(R+1,A+1,2*p1),1e-300,1)
        logV[:,:,states['F']] = (D>0.5)*( numpy.log(safe_fixed_beta) - numpy.log(2*p1) ) 
        
        # Polymorphic 
        logV[:,:,states['P']] = (D>0.5)*( 0.0 )
           
        logVmax = logV.max(axis=2)
        
        for l in xrange(0,L):
            logV[:,:,l] -= logVmax
                
        V = numpy.exp(logV)
        
        # Compute forward and backward tables using recursion relations in SI
        
        # All mutations start in ancestral state
        F[:,0,states['A']] = 1.0
        F[:,0,states['E']] = 0.0
        F[:,0,states['F']] = 0.0
        F[:,0,states['P']] = 0.0
        
        # Mutations cannot end in ancestral state
        B[:,T-1,states['A']]=0.0
        B[:,T-1,states['E']] = 1.0/(L-1)
        B[:,T-1,states['F']] = 1.0/(L-1)
        B[:,T-1,states['P']] = 1.0/(L-1)
        
        
        for t in xrange(1,T):
            F[:,t,:] = V[:,t,:]*( numpy.dot( F[:,t-1,:], numpy.transpose(Q[t-1,:,:])) )
            sum_F = F[:,t,:].sum(axis=1)
            sum_F = sum_F + (sum_F < 1e-12)
            F[:,t,:] /= numpy.outer(sum_F, numpy.ones(L))    
            
        for t in reversed(xrange(0,T-1)):
            B[:,t,:] = numpy.dot(B[:,t+1,:]*V[:,t+1,:], Q[t,:,:])
            sum_B = B[:,t,:].sum(axis=1)
            sum_B = sum_B + (sum_B < 1e-12)
            B[:,t,:] /= numpy.outer(sum_B, numpy.ones(L))    
          
        # calculate probability of being in a state
        # and probability of transition
        
        # probability of state
        for t in xrange(0,T):
            Pstate[:,t,:] = F[:,t,:]*B[:,t,:]
            sum_P = Pstate[:,t,:].sum(axis=1)
            sum_P = sum_P+(sum_P<1e-12)
            Pstate[:,t,:] /= numpy.outer(sum_P, numpy.ones(L))    
        
        ancestral_error_numerators = (A*Pstate[:,:,0]).sum(axis=1)+(A*Pstate[:,:,1]).sum(axis=1)
        ancestral_error_denominators = (D*Pstate[:,:,0]).sum(axis=1)+(D*Pstate[:,:,1]).sum(axis=1)
        derived_error_numerators = (R*Pstate[:,:,2]).sum(axis=1)
        derived_error_denominators = (D*Pstate[:,:,2]).sum(axis=1)
        
        p0 = ancestral_error_numerators/(ancestral_error_denominators+(ancestral_error_denominators<1e-12))
        p0 = numpy.clip(p0,min_error,epsE)
        p0 = numpy.outer(p0,numpy.ones(T))
        
            
        p1 = derived_error_numerators/(derived_error_denominators+(derived_error_denominators<1e-12))
        p1 = numpy.clip(p1,min_error,epsF)
        p1 = numpy.outer(p1,numpy.ones(T))
        
    
        if current_iteration==num_iterations:
            # infer most likely sequence of states
    
            max_likelihood_table[:,0,:] = F[:,0,:]
            
            for t in xrange(1,T):
            
                new_likelihoods = max_likelihood_table[:,t-1,None,:]*Q[None,t-1,:,:]
                
                max_likelihood_table[:,t,:] = new_likelihoods.max(axis=2)
                state_pointer_table[:,t-1,:] = new_likelihoods.argmax(axis=2)
                
                if t==T-1:
                    # not possible to be in unborn state at end
                    max_likelihood_table[:,t,:] *= B[:,T-1,:]
                
                max_likelihood_table[:,t,:] *= V[:,t,:]                
                                
                scale = max_likelihood_table[:,t,:].max(axis=1)
                scale = numpy.outer(numpy.clip(scale,1e-300,1e300),numpy.ones(L))
                max_likelihood_table[:,t,:] /= scale
                
            # follow path back
            Ls = numpy.zeros((M,T),dtype=numpy.int32)
            
            Ls[:,-1] = max_likelihood_table[:,-1,:].argmax(axis=1)
            for t in reversed(xrange(1,T)):
                for m in xrange(0,M):
                    Ls[m,t-1] = state_pointer_table[m,t-1,Ls[m,t]]
        
    return Pstate,Ls,p0,p1

