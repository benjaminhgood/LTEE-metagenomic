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

states = parse_file.clade_hmm_states

# to have no clades, set f0s to all zeros

def infer_hmm(A, D, f0s=None, num_iterations=10, penalize_clades=True, infer_fs=True, debug=False):
    # run baum welch algorithm
    
    if f0s==None:
        estimate_f0 = True
    else:
        estimate_f0 = False
        fmajor0s, fminor0s = f0s
    
    # fixed params
    epsF = 0.025 # error rate from 100% for fixed mutations
    epsE = 0.025 # error rate from 0% for extinct or unborn mutations
    min_error = 1e-03 # illumina error rate
    eps = 0.025
    min_clade_freq = min_error # was 1e-02
    
    # transition rates for markov chain    
    overall_birth_rate = 2e-02
    ancestral_birth_rate = overall_birth_rate/4 
    clade_birth_rate = overall_birth_rate*3/8
    ancestral_fix_rate = 0.5
    clade_fix_rate = 0.5
    rebirth_rate = 1e-06
    clade_jump_rate = 1e-20
        
    M = A.shape[0] # mutations
    T = A.shape[1]     # timepoints
    ts = numpy.arange(0,T)
    
    L = len(states.keys())
    
    fmajors = numpy.zeros_like(ts)
    fminors = numpy.zeros_like(ts)
    
    R = D-A
    safe_D = D+(D<0.5) 
    safe_f = A/safe_D
    safe_f = numpy.clip(safe_f, min_clade_freq, 1-min_clade_freq)
    
    p0 = numpy.ones_like(safe_f)*1e-02
    p1 = numpy.ones_like(safe_f)*1e-02
    
    if estimate_f0:
        zeroth_fs = (A/(D+(D==0)))
        # calculate a measure of "heterozygosity" as a function of time
        Hs = ((zeroth_fs>0.5)*(zeroth_fs<0.8)*((zeroth_fs-0.5)*(0.8-zeroth_fs))).sum(axis=0)
        #Hs = ((zeroth_fs>0.2)*(zeroth_fs<0.8)*(zeroth_fs-0.2)*(0.8-zeroth_fs)).sum(axis=0)
        # find point where this is maximized
        tstar = Hs.argmax()
        
        sys.stderr.write("Diversity peak at %g\n" % tstar)
        
        # Get indices that belong to each cluster
        # This may include *some* misclassification, but hopefully a minority 
        # so that "typical" properties are the same
        putative_majority_idxs = (zeroth_fs[:,tstar]<0.8)*(zeroth_fs[:,tstar]>0.5)
        putative_minority_idxs = (zeroth_fs[:,tstar]>0.2)*(zeroth_fs[:,tstar]<0.5)
        
         
        fmajors = numpy.zeros_like(Hs)
        fminors = numpy.zeros_like(Hs)
        good_majority_idxs = numpy.zeros(T,dtype=numpy.bool_)
        good_minority_idxs = numpy.zeros(T,dtype=numpy.bool_)
        
        # for timepoints beyond tstar, use the inferred collection of (fixed) mutations
        # to estimate the clade frequencies
        for t in xrange(tstar,T):
            
            fs = zeroth_fs[putative_majority_idxs*(D[:,t]>0.5),t]
            if len(fs) < 1:
                pass
            else:
                fmajors[t]=numpy.median(fs)
                good_majority_idxs[t]=True
                fminors[t] = 1-fmajors[t]
                good_minority_idxs[t]=True
            
            
            #fs = zeroth_fs[putative_minority_idxs*(D[:,t]>0.5),t]
            #if len(fs) < 1:
            #    pass
            #else:
            #    fminors[t]=numpy.median(fs)
            #    good_minority_idxs[t]=True
        
        
        # if the majority clade at time tmax becomes the minority at the end, switch them
        if fmajors[good_majority_idxs][-1] < fminors[good_minority_idxs][-1]:
            putative_majority_idxs, putative_minority_idxs = putative_minority_idxs, putative_majority_idxs
            good_majority_idxs, good_minority_idxs = good_minority_idxs, good_majority_idxs
            fmajors,fminors = fminors,fmajors
        
        # now that we have an estimate of the clade trajectories after tmax, 
        # we can use these to re-classify the mutations at tmax into the proper clades
        for i in xrange(0,len(putative_majority_idxs)):
            
            if putative_majority_idxs[i]>0 or putative_minority_idxs[i]>0:
                
                D1 = numpy.square(zeroth_fs[i,(D[i,:]>0)*(good_majority_idxs)]-fmajors[(D[i,:]>0)*(good_majority_idxs)] ).mean()
                D2 = numpy.square(zeroth_fs[i,(D[i,:]>0)*(good_minority_idxs)]-fminors[(D[i,:]>0)*(good_minority_idxs)] ).mean()
                D3 = numpy.square(zeroth_fs[i,(D[i,:]>0)*(ts>=tstar)]-0).mean()
                D4 = numpy.square(zeroth_fs[i,(D[i,:]>0)*(ts>=tstar)]-1).mean()

                if D1 < min([D2,D3,D4]):
                    # actually a major
                    putative_majority_idxs[i] = True
                    putative_minority_idxs[i] = False
                else:
                    # not really in the major clade
                    putative_majority_idxs[i]=False
                    
                    if D2 < min([D1,D3,D4]):
                        # actually in the minor clade
                        putative_minority_idxs[i]=True
                    else:
                        putative_minority_idxs[i]=False
        
        # now we have a decent putative majority set
        # so we can use these to estimate the clade frequency at all times
        for t in xrange(0,T):
            fs = zeroth_fs[putative_majority_idxs*(D[:,t]>0.5), t]
            if len(fs) > 0:
                fmajors[t] = fs[fs>=(fs.max()/2.0)].mean() 
                good_majority_idxs[t]=True
            else:
                good_majority_idxs[t]=False  
                  
            fs = zeroth_fs[putative_minority_idxs*(D[:,t]>0.5), t]
            if len(fs) > 0:
                fminors[t] = fs[fs>=(fs.max()/2.0)].mean()
                good_minority_idxs[t]=True
            else:
                good_minority_idxs[t]=False
                
            if fmajors[t]+fminors[t]>1:
                fminors[t]=1-fmajors[t]
        
        # do interpolation to fill in missing timepoints
        interpolated_fmajors = interp1d(ts[good_majority_idxs],fmajors[good_majority_idxs], bounds_error=False,fill_value=1e-02)
        fmajors = interpolated_fmajors(ts)
        interpolated_fminors = interp1d(ts[good_minority_idxs],fminors[good_minority_idxs], bounds_error=False,fill_value=1e-02)
        fminors = interpolated_fminors(ts)
        #fminors = numpy.fmin(fminors,1-fmajors)
        
        #print ", ".join([str(f) for f in fmajors])
        #print ", ".join([str(f) for f in fminors])
        
    else:
        fmajors = fmajor0s
        fminors = fminor0s
    
    if ((fmajors+fminors)>0.2).sum() < 1:
        print "Disallowing transitions!"
        # do not allow transitions to clades
        #print "disallowing transitions"
        ancestral_birth_rate = overall_birth_rate 
        clade_birth_rate = 0
    
    
    ftotals = numpy.clip(fmajors+fminors,min_clade_freq,1-min_clade_freq)    
    fmajors = numpy.clip(fmajors,min_clade_freq,1-min_clade_freq)
    fminors = numpy.clip(fminors,min_clade_freq,1-min_clade_freq)
    
    #print ", ".join(["%0.2f" % f for f in fmajors])
    
    Q = numpy.zeros((T,L,L))
    F = numpy.zeros((M,T,L))*1.0
    B = numpy.zeros((M,T,L))*1.0
    
    Pstate = numpy.zeros((M,T,L))*1.0
    Ptrans = numpy.zeros((M,T,L,L))*1.0
    
    max_likelihood_table = numpy.zeros((M,T,L))*1.0
    state_pointer_table = numpy.zeros((M,T,L),dtype=numpy.int32)
    
        
    for current_iteration in xrange(1,num_iterations+1):
    
        #print fmajors[-1], fminors[-1]
    
        sys.stderr.write("Iteration %d\n" % current_iteration)
        
        # calculate Qs and Vs
        #sys.stderr.write("Calculating Q,V...\n")
        
        Q = numpy.zeros((T, L, L))
        
        for t in xrange(0,T):
            
            # Ancestral
            #     to          from
            Q[t,states['A'],states['A']] = 1-overall_birth_rate
            Q[t,states['E'],states['E']] = 1-rebirth_rate
            Q[t,states['A'],states['E']] = rebirth_rate
            
            # Extinct
            Q[t,states['E'],states['PB']] = ancestral_fix_rate/2
            Q[t,states['E'],states['PM']] = clade_fix_rate/2
            Q[t,states['E'],states['Pm']] = clade_fix_rate/2
            
            # Fixed basal
            Q[t,states['FB'],states['FB']] = 1
            Q[t,states['FB'],states['PB']] = ancestral_fix_rate/2
            Q[t,states['FB'],states['PB*']] = ancestral_fix_rate
            
            # Fixed minority
            Q[t,states['Fm'],states['Fm']] = 1-clade_jump_rate
            Q[t,states['Fm'],states['Pm']] = clade_fix_rate/2
            
            # Fixed majority 
            Q[t,states['FM'],states['FM']] = 1-clade_jump_rate
            Q[t,states['FM'],states['PM']] = clade_fix_rate/2
            
            # Polymorphic basal
            Q[t,states['PB'],states['A']] = ancestral_birth_rate
            #Q[t,states['PB'],states['E']] = rebirth_rate/3
            Q[t,states['PB'],states['PB']] = 1 - ancestral_fix_rate
            
            # Poly-minority
            Q[t,states['Pm'],states['A']] = clade_birth_rate
            #Q[t,states['Pm'],states['E']] = rebirth_rate/3
            Q[t,states['Pm'],states['Pm']] = 1-clade_fix_rate
            
            # Poly-majority
            Q[t,states['PM'],states['A']] = clade_birth_rate
            #Q[t,states['PM'],states['E']] = rebirth_rate/3
            Q[t,states['PM'],states['PM']] = 1-clade_fix_rate
            
            # From within-clade to population-wide fixation
            Q[t,states['PB*'],states['PB*']] = 1-ancestral_fix_rate
            Q[t,states['PB*'],states['FM']] = clade_jump_rate
            Q[t,states['PB*'],states['Fm']] = clade_jump_rate
            
        
        logV = numpy.zeros((M,T,L))*1.0
        V = numpy.zeros_like(logV)
        
        #sys.stderr.write("Calculating emission probabilities...\n")
        
        bigfmajors = numpy.outer(numpy.ones(M),fmajors)
        bigfminors = numpy.outer(numpy.ones(M),fminors)
        bigftotals = numpy.outer(numpy.ones(M),ftotals)
    
        #bigfmajors += (p0-bigfmajors)*(bigfmajors<p0)
        #bigfmajors += (1-p1-bigfmajors)*(bigfmajors>(1-p1))
        
        #bigfminors += (p0-bigfminors)*(bigfminors<p0)
        #bigfminors += (1-p1-bigfminors)*(bigfminors>(1-p1))
        
        #bigftotals += (p0-bigftotals)*(bigftotals<p0)
        #bigftotals += (1-p1-bigftotals)*(bigftotals>(1-p1))
           
        
        #fpfixeds_lower_error = numpy.clip(bigftotals-eps,1e-03,1-1e-03)
        #fpminors_upper_error = numpy.clip(bigfminors+eps,1e-03,1-1e-03)
        #fpminors_lower_error = numpy.clip(bigfminors-eps,1e-03,1-1e-03)
        #fpmajors_upper_error = numpy.clip(bigfmajors+eps,1e-03,1-1e-03)
        #fpmajors_lower_error = numpy.clip(bigfmajors-eps,1e-03,1-1e-03)
        
        fpfixeds_lower_error = numpy.clip(bigftotals-p0,1e-03,1-1e-03)
        fpminors_upper_error = numpy.clip(bigfminors+p0,1e-03,1-1e-03)
        fpminors_lower_error = numpy.clip(bigfminors-p0,0, 1-1e-03)
        fpmajors_upper_error = numpy.clip(bigfmajors+p0,1e-03,1-1e-03)
        fpmajors_lower_error = numpy.clip(bigfmajors-p0,1e-03,1-1e-03)
        
        
        # Cap at error thresholds
        fpmajors_lower_error += (0-fpmajors_lower_error)*(fpmajors_upper_error<(2*p0))
        fpmajors_upper_error += (2*p0 - fpmajors_upper_error)*(fpmajors_upper_error<(2*p0))
        
        fpminors_lower_error += (0-fpminors_lower_error)*(fpminors_upper_error<(2*p0))
        fpminors_upper_error += (2*p0 - fpminors_upper_error)*(fpminors_upper_error<(2*p0))
          
        fpfixeds_lower_error += (1-2*p1-fpfixeds_lower_error)*(fpfixeds_lower_error>(1-2*p1))
            
        if penalize_clades:
            #fixed_ancestral_lower_bound = bigftotals
            #polymorphic_major_upper_bound = bigfmajors
            #polymorphic_minor_upper_bound = bigfminors    
            fixed_ancestral_lower_bound = fpfixeds_lower_error
            polymorphic_major_upper_bound = fpmajors_upper_error
            polymorphic_minor_upper_bound = fpminors_upper_error
        
        else:
            fixed_ancestral_lower_bound = numpy.zeros_like(bigftotals)
            polymorphic_major_upper_bound = numpy.ones_like(bigfmajors)
            polymorphic_minor_upper_bound = numpy.ones_like(bigfminors)  
        
        
        # Unborn   
        safe_unborn_beta = numpy.clip(incomplete_beta(A+1,R+1,2*p0),1e-300,1)
        logV[:,:,states['A']] = (D>0.5)*( numpy.log(safe_unborn_beta) - numpy.log(2*p0) ) 
        
        # Extinct 
        safe_extinct_beta = numpy.clip(incomplete_beta(A+1,R+1,2*p0),1e-300,1)
        logV[:,:,states['E']] = (D>0.5)*( numpy.log(safe_extinct_beta) - numpy.log(2*p0) )
        
        # fixed ancestral
        safe_fixed_beta = numpy.clip(incomplete_beta(R+1,A+1,2*p1),1e-300,1)
        logV[:,:,states['FB']] = (D>0.5)*( numpy.log(safe_fixed_beta) - numpy.log(2*p1) )
        
        # fixed minority
        safe_fixed_minor_beta = numpy.clip(incomplete_beta(A+1,R+1,fpminors_upper_error)-incomplete_beta(A+1,R+1,fpminors_lower_error),1e-300,1)
        logV[:,:,states['Fm']] = (D>0.5)*( numpy.log(safe_fixed_minor_beta) - numpy.log(fpminors_upper_error-fpminors_lower_error) )
        
        # fixed majority
        safe_fixed_major_beta = numpy.clip(incomplete_beta(A+1,R+1,fpmajors_upper_error)-incomplete_beta(A+1,R+1,fpmajors_lower_error),1e-300,1)
        logV[:,:,states['FM']] = (D>0.5)*( numpy.log(safe_fixed_major_beta) - numpy.log(fpmajors_upper_error-fpmajors_lower_error) ) 
        
        # polymorphic ancestral 
        safe_polymorphic_ancestral_beta = numpy.clip(incomplete_beta(R+1,A+1, 1-fixed_ancestral_lower_bound),1e-300,1)
        logV[:,:,states['PB']] = (D>0.5)*( numpy.log(safe_polymorphic_ancestral_beta) - numpy.log(1-fixed_ancestral_lower_bound) )
        
        # poly minority
        safe_polymorphic_minor_beta = numpy.clip(incomplete_beta(A+1,R+1, polymorphic_minor_upper_bound),1e-300,1)
        logV[:,:,states['Pm']] = (D>0.5)*( numpy.log(safe_polymorphic_minor_beta) - numpy.log(polymorphic_minor_upper_bound) )
        
        # poly majority
        safe_polymorphic_major_beta = numpy.clip(incomplete_beta(A+1,R+1, polymorphic_major_upper_bound),1e-300,1)
        logV[:,:,states['PM']] = (D>0.5)*( numpy.log(safe_polymorphic_major_beta) - numpy.log(polymorphic_major_upper_bound) )
        
        # polymorphic ancestral but then fixed 
        safe_polymorphic_ancestral_beta = numpy.clip(incomplete_beta(R+1,A+1, 1-fixed_ancestral_lower_bound),1e-300,1)
        logV[:,:,states['PB*']] = (D>0.5)*( numpy.log(safe_polymorphic_ancestral_beta) - numpy.log(1-fixed_ancestral_lower_bound) )
        
          
        logVmax = logV.max(axis=2)
        
        for l in xrange(0,L):
            logV[:,:,l] -= logVmax
                
        V = numpy.exp(logV)
        
        #if debug==True:
        #    for state in states.keys():
        #        print state, "V", ", ".join(["%0.2e" % v for v in V[0,-25:-15,states[state]]])
            
        #    print ", ".join(["%0.3f" % f for f in safe_extinct_beta[0,-25:-15]])
        #    print ", ".join(["%0.3f" % f for f in safe_fixed_minor_beta[0,-25:-15]])
            
        #    print ", ".join(["%0.3f" % f for f in (fpminors_upper_error-fpminors_lower_error)[0,-25:-15]])
        #    print ", ".join(["%0.3f" % (2*f) for f in p0[0,-25:-15]])
        
        # Compute forward and backward tables
        
        # Always start in unborn state
        F[:,0,states['A']] = 1.0 
        F[:,0,1:] = 0.0 
        
        # Must finish in some kind of fixed state
        
        B[:,T-1,states['A']] = 0.0
        B[:,T-1,states['Pm']] = 0.0
        B[:,T-1,states['PB*']] = 0.0
        
        Lf = 6.0
        B[:,T-1,states['PB']] = 1/Lf
        B[:,T-1,states['PM']] = 1/Lf
        B[:,T-1,states['E']] = 1/Lf
        B[:,T-1,states['FB']] = 1/Lf
        B[:,T-1,states['FM']] = 1/Lf
        B[:,T-1,states['Fm']] = 1/Lf
              
        for t in xrange(1,T):
            F[:,t,:] = V[:,t,:]*( numpy.dot( F[:,t-1,:], numpy.transpose(Q[t-1,:,:])) )
            max_F = numpy.outer(F[:,t,:].max(axis=1), numpy.ones(L))
            sum_F = (F[:,t,:]/max_F).sum(axis=1)
            F[:,t,:] /= numpy.outer(sum_F, numpy.ones(L))*max_F   
            
        for t in reversed(xrange(0,T-1)):
            B[:,t,:] = numpy.dot(B[:,t+1,:]*V[:,t+1,:], Q[t,:,:])
            max_B = numpy.outer(B[:,t,:].max(axis=1), numpy.ones(L))
            sum_B = (B[:,t,:]/max_B).sum(axis=1)
            B[:,t,:] /= numpy.outer(sum_B, numpy.ones(L))*max_B
        
        F = numpy.clip(F,1e-300,1)
        B = numpy.clip(B,1e-300,1)
          
        # calculate probability of being in a state
        # and probability of transition
        #sys.stderr.write("Calculating state probabilities...\n")
        # probability of state
        for t in xrange(0,T):
            Pstate[:,t,:] = F[:,t,:]*B[:,t,:]
            max_P = numpy.outer(Pstate[:,t,:].max(axis=1), numpy.ones(L))
            max_P = numpy.clip(max_P,1e-300,1)
            sum_P = (Pstate[:,t,:]/max_P).sum(axis=1)
            Pstate[:,t,:] /= numpy.outer(sum_P, numpy.ones(L)) * max_P   
        
        #if debug==True:
        #    for state in states.keys():
        #        print state, "B", ", ".join(["%0.2e" % p for p in B[0,-25:-15,states[state]]])
        #        print state, "F", ", ".join(["%0.2e" % p for p in F[0,-25:-15,states[state]]])
        
        #print p0[0,-1]
        #print p1[0,-1]
        #print Pstate[0,-1,:]
                       
        # update fs and deltas
        
        #print "Iteration", current_iteration
        #print Pstate
        
        # infer p0 and p1
        
        ancestral_error_numerators = (A*Pstate[:,:,states['A']]).sum(axis=1)+(A*Pstate[:,:,states['E']]).sum(axis=1)
        ancestral_error_denominators = (D*Pstate[:,:,states['A']]).sum(axis=1)+(D*Pstate[:,:,states['A']]).sum(axis=1)
        derived_error_numerators = (R*Pstate[:,:,states['FB']]).sum(axis=1)
        derived_error_denominators = (D*Pstate[:,:,states['FB']]).sum(axis=1)
        
        p0 = ancestral_error_numerators/(ancestral_error_denominators+(ancestral_error_denominators<1e-12))
        p0 = numpy.clip(p0,min_error,epsE)
        p0 = numpy.outer(p0,numpy.ones(T))
            
        p1 = derived_error_numerators/(derived_error_denominators+(derived_error_denominators<1e-12))
        p1 = numpy.clip(p1,min_error,epsF)
        p1 = numpy.outer(p1,numpy.ones(T))
        
        if infer_fs:
        
            # fmajors
            # decode ML states
            #Ls = numpy.zeros_like(A) 
            #for i in xrange(0,A.shape[0]):
            #    Ls[i] = (Pstate[i,:,:]).argmax(axis=1)
            
            major_poly_numerators = numpy.zeros(T)*1.0
            major_poly_denominators = numpy.zeros_like(major_poly_numerators)
            minor_poly_numerators = numpy.zeros_like(major_poly_numerators)
            minor_poly_denominators = numpy.zeros_like(major_poly_numerators)
            
            final_Ls = Pstate[:,-1,:].argmax(axis=1)
            
            majority_idxs = (final_Ls==states['FM'])
            minority_idxs = (final_Ls==states['Fm'])
            
            major_numerators = (A*Pstate[:,:,states['FM']]*majority_idxs[:,None]).sum(axis=0)+major_poly_numerators
            major_denominators = (D*Pstate[:,:,states['FM']]*majority_idxs[:,None]).sum(axis=0)+major_poly_denominators
            minor_numerators = (A*Pstate[:,:,states['Fm']]*minority_idxs[:,None]).sum(axis=0)+minor_poly_numerators
            minor_denominators = (D*Pstate[:,:,states['Fm']]*minority_idxs[:,None]).sum(axis=0)+minor_poly_denominators
            
            fmajors = major_numerators/(major_denominators+(major_denominators<1e-12))
            fminors = minor_numerators/(minor_denominators+(minor_denominators<1e-12))
            
            bad_indices = ((fmajors+fminors)>1)
        
            fminors[bad_indices] = 1-fmajors[bad_indices]
        
            fmajors = numpy.clip(fmajors, min_clade_freq, 1-min_clade_freq)
            fminors = numpy.clip(fminors, min_clade_freq, 1-min_clade_freq)
        
            if fmajors[-1] < fminors[-1]:
                print "Flipping!"
                fmajors, fminors = fminors, fmajors

            #print fs
            
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
        
        
    
    if infer_fs:
    
        # Finishing touches to interpolate fmajor and fminor
        
        final_Ls = Pstate[:,-1,:].argmax(axis=1)
     
        majority_idxs = (final_Ls==states['FM'])
        minority_idxs = (final_Ls==states['Fm'])
        nmajors = (Pstate[majority_idxs,:,states['FM']]+Pstate[majority_idxs,:,states['PM']]).sum()   
        nminors = (Pstate[majority_idxs,:,states['FM']]+Pstate[majority_idxs,:,states['PM']]).sum()   
     
        good_idxs = major_denominators>1
        good_idxs[0] = 1
        good_idxs[-1] = 1
        # interpolate fs
        interpolated_fs = interp1d(ts[good_idxs],fmajors[good_idxs])
        fmajors = interpolated_fs(ts)
        fmajors = fmajors*(nmajors>0.49)
    
        good_idxs = minor_denominators>1
        good_idxs[0] = 1
        good_idxs[-1] = 1
        # interpolate fs
        interpolated_fs = interp1d(ts[good_idxs],fminors[good_idxs])
        fminors = interpolated_fs(ts)
        fminors = fminors*(nminors>0.49)

    return Pstate,Ls,fmajors,fminors

