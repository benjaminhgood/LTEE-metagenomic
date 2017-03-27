import numpy
from scipy.interpolate import interp1d
import parse_file

def calculate_loglikelihood_ratio(A0,D0,A1,D1):
    
    R0 = D0-A0
    R1 = D1-A1
    
    logL = 0
    logL += (A1*numpy.log((A0+A1+(A1==0))/(D0+D1+(A1==0))*(D1+(A1==0))/(A1+(A1==0)))).sum()
    logL += (R1*numpy.log((R0+R1+(R1==0))/(D0+D1+(R1==0))*(D1+(R1==0))/(R1+(R1==0)))).sum()
    logL += (A0*numpy.log((A0+A1+(A0==0))/(D0+D1+(A0==0))*(D0+(A0==0))/(A0+(A0==0)))).sum() 
    logL +=         (R0*numpy.log((R0+R1+(R0==0))/(D0+D1+(R0==0))*(D0+(R0==0))/(R0+(R0==0)))).sum()

    dof = ((D0>0.5)*(D1>0.5)).sum()*1.0
    
    if dof<0.5:
        return -1,dof
    else:
        return -2*logL/dof, dof
    

###########################
#
# Estimates the slope of the x vector
# using linear regression
#
###########################
def estimate_slope(ts, xs):
        return ((xs*ts).mean()-xs.mean()*ts.mean())/((ts*ts).mean()-ts.mean()*ts.mean())
    

############################
#
# Takes a list of (ts,xs) pairs
# (where ts may partially overlap)
# and returns estimate of avg xs 
# at timepoints where this can be estimated
#
############################    
def average_trajectories(trajectories):

    avg_map = {}
    
    for ts,xs in trajectories:
        for t,x in zip(ts,xs):
            if t not in avg_map:
                avg_map[t] = {'x': 0.0, 'n': 0.0}
            avg_map[t]['x']+=x
            avg_map[t]['n']+=1
    
    for t in avg_map.keys():
        avg_map[t]['x'] /= avg_map[t]['n']
    
    avg_ts = []
    avg_xs = []
    for t in sorted(avg_map.keys()):
        avg_ts.append(t)
        avg_xs.append(avg_map[t]['x'])
        
    return numpy.array(avg_ts), numpy.array(avg_xs)
    

    
########
#
# Takes a trajectory and zeros out alt and depth entries where coverage
# is less than minimum threshold, or deletion is inferred to have happened
#
########
def mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue, min_depth=parse_file.default_min_depth):
    
    # first make a copy of alts and depths 
    # so that we can modify in place without worrying
    masked_alts = numpy.copy(alts)
    masked_depths = numpy.copy(depths)
    
    # zero out timepoints that don't pass depth threshold
    masked_alts[masked_depths<min_depth] = 0
    masked_depths[masked_depths<min_depth] = 0
    
    #masked_alts -= masked_alts*(masked_depths < min_depth)
    #masked_depths -= masked_depths*(masked_depths < min_depth)
    
    # did we infer that a deletion happened?
    #if (var_type=='sv' and depth_fold_change < -2 and depth_change_pvalue < 1e-04) or (var_type!='sv' and depth_fold_change < -1 and depth_change_pvalue < 1e-03):
    if depth_change_pvalue < 1e-02:   
       # deletion nearby, trim timecourse
       masked_alts[cutoff_idx:] = 0
       masked_depths[cutoff_idx:] = 0
       
    good_idxs = numpy.nonzero(masked_depths>0.5)[0]
    
    return good_idxs, masked_alts, masked_depths

###########
#
# Naive frequency estimator, # alts / # depths (0 if no depth)
#
###########
def estimate_frequencies(alts,depths):
    return alts*1.0/(depths+(depths==0))

###########
#
# Naive frequency estimator for clones, restricted to clones with >=min_depth
#
###########
def estimate_clone_frequencies(times, alts, depths, min_depth=20, allowed_times=None):

    masked_times = times[depths>min_depth]
    masked_alts = alts[depths>min_depth]
    masked_depths = depths[depths>min_depth]
    masked_freqs = masked_alts*1.0/masked_depths
    
    return masked_times, masked_freqs
    
    
###########
#
# Estimate depth fold change (log2) relative to genome-wide median
#
###########    
def estimate_depth_fold_changes(avg_depths, depths, min_depth=20):
    
    if (depths<=0).all():
        return numpy.zeros_like(depths)
    
    normalization = depths[depths>min_depth][0]/avg_depths[depths>min_depth][0]
    
    return numpy.log2((depths+(depths==0))/(avg_depths+(avg_depths==0))/normalization)

###########
#
# Creates a continuous interpolation function for frequency trajectory,
# so that it can be evaluated anywhere on time interval
#
###########    
def create_interpolation_function(times,freqs,tmax=100000,kind='linear'):
    # can create it for anything!

    padded_times = numpy.zeros(len(times)+1)
    padded_freqs = numpy.zeros(len(times)+1)
    padded_times[0:len(times)] = times
    padded_freqs[0:len(times)] = freqs
    padded_times[-1] = tmax
    padded_freqs[-1] = freqs[-1]
    
    #xis = numpy.log((alts+1e-01)/(depths-alts+1e-01))
    #interpolating_function = interp1d(times, xis, kind='linear',bounds_error=False)
    #interpolated_xis = interpolating_function(theory_times)
    #interpolated_freqs = 1.0/(1+numpy.exp(-interpolated_xis))
    #freqs = 1.0/(1+numpy.exp(-xis))
    
    interpolating_function = interp1d(padded_times, padded_freqs, kind=kind,bounds_error=True)
    
    return interpolating_function

def calculate_appearance_time_from_hmm(times,fs,Ls):
    
    # Calculate time tstar at which f(t) is largest
    if Ls[-1]==parse_file.FIXED:
        # If fixed, this is final timepoint
        tstar = times[-1]
    else:
        # Otherwise, pick point where f(t) is largest
        # (restrict to polymorphic timepoints so that we
        #  don't focus on an error fluctuation)
        tstar = (times[Ls==parse_file.POLYMORPHIC])[fs[Ls==parse_file.POLYMORPHIC].argmax()]

    #first_polymorphic_idx = numpy.nonzero((Ls==parse_file.UNBORN)*(times<tstar))[0][-1]+1
    first_polymorphic_idx = numpy.nonzero( numpy.logical_or(Ls==parse_file.UNBORN, Ls==parse_file.EXTINCT)*(times<tstar) )[0][-1]+1
    
    T0 = times[first_polymorphic_idx]-250
    
    return T0

def calculate_appearance_fixation_time_from_hmm(times,fs,Ls):
    
    extinct_idxs = numpy.array([l in parse_file.well_mixed_extinct_states for l in Ls])
    polymorphic_idxs = numpy.array([l in parse_file.well_mixed_polymorphic_states for l in Ls])
    non_polymorphic_idxs = numpy.logical_not(polymorphic_idxs)
    num_polymorphic_idxs = (polymorphic_idxs).sum()
    
    # Calculate time tstar at which f(t) is largest
    if Ls[-1] in parse_file.well_mixed_fixed_states:

        # If fixed, this is final timepoint
        tstar = times[-1]

    else:
        # Otherwise, pick point where f(t) is largest
        # (restrict to polymorphic timepoints so that we
        #  don't focus on an error fluctuation)
        
        tstar = (times[polymorphic_idxs])[fs[polymorphic_idxs].argmax()]

       
    appearance_time = times[numpy.nonzero((times<=tstar)*extinct_idxs)[0][-1]]+250
    
    later_non_polymorphic_idxs = numpy.nonzero((times>=appearance_time)*non_polymorphic_idxs)[0]
    
    
    if len(later_non_polymorphic_idxs) == 0:
        # Stays polymorphic until the end
        fixation_time = 1000000 # never fixes
        transit_time = times[-1]+250 - appearance_time
    else:
        fixation_time = times[later_non_polymorphic_idxs[0]]-250
        transit_time = fixation_time - appearance_time
        
    return appearance_time, fixation_time, transit_time

def calculate_appearance_fixation_time_from_clade_hmm(times,fs,Ls):
    
    extinct_idxs = numpy.array([l in parse_file.clade_extinct_states for l in Ls])
    polymorphic_idxs = numpy.array([l in parse_file.clade_polymorphic_states for l in Ls])
    non_polymorphic_idxs = numpy.logical_not(polymorphic_idxs)
    num_polymorphic_idxs = (polymorphic_idxs).sum()
    
    # Calculate time tstar at which f(t) is largest
    if Ls[-1] in parse_file.clade_fixed_states:

        # If fixed, this is final timepoint
        tstar = times[-1]

    else:
        # Otherwise, pick point where f(t) is largest
        # (restrict to polymorphic timepoints so that we
        #  don't focus on an error fluctuation)
        
        tstar = (times[polymorphic_idxs])[fs[polymorphic_idxs].argmax()]

       
    appearance_time = times[numpy.nonzero((times<=tstar)*extinct_idxs)[0][-1]]+250
    
    later_non_polymorphic_idxs = numpy.nonzero((times>=appearance_time)*non_polymorphic_idxs)[0]
    
    
    if len(later_non_polymorphic_idxs) == 0:
        # Stays polymorphic until the end
        fixation_time = 1000000 # never fixes
        transit_time = times[-1]+250 - appearance_time
    else:
        fixation_time = times[later_non_polymorphic_idxs[0]]-250
        transit_time = fixation_time - appearance_time
    
    return appearance_time, fixation_time, transit_time

def split_clade_hmm(times,fs,Ls):
    # splits into independent mutations
    
    independent_runs = []
    
    previous_state = parse_file.clade_hmm_states['A']
    in_run = False
    current_run = []
    for t_idx in xrange(0,len(times)):
        current_state = Ls[t_idx]
        if (not in_run) and previous_state==parse_file.clade_hmm_states['A'] and current_state!=parse_file.clade_hmm_states['A']:
            # transitioned from ancestral to polymorphic!
            # start a new run!
            #print "Starting new run!", t_idx
            in_run = True
            
        if in_run:
            current_run.append(t_idx)
            
        if in_run and (current_state in parse_file.clade_fixed_states or current_state==parse_file.clade_hmm_states['E'] or current_state==parse_file.clade_hmm_states['A']):
            # The mutation has fixed or gone extinct, so we are done with the run!
            in_run = False
            #print current_state, previous_state, in_run, current_run
            independent_runs.append( numpy.array(current_run) )
            current_run = []
            
        previous_state = current_state
        
    if len(current_run) > 0:
        # we were in the middle of a run when experiment ended
        independent_runs.append( numpy.array(current_run) )
       
    return independent_runs
    
        
def calculate_appearance_time_from_freqs(times, freqs):
    
    # two ways to do it
    # first time above 10%:
    if (freqs>=0.1).any():
        
        # first time above fmin
        #t = times[numpy.nonzero(freqs>=0.1)[0][0]]
        
        # or last time below fmin 
        # we've been using this one
        max_idx = freqs.argmax()    
        zero_idxs = numpy.nonzero(freqs[:max_idx]<=0.1)[0]
            
        if len(zero_idxs) == 0:
            t = times[max_idx]
        else:
            t = times[zero_idxs[-1]+1]
            
    else:
        t = times[freqs.argmax()]
    
    return t

def calculate_appearance_time(ts, fs, state_Ls, Ls):
    # Three choices: either from freqs, simple HMM, or clade HMM
    #return calculate_appearance_time_from_freqs(ts,fs)
    #return calculate_appearance_time_from_hmm(ts,fs,state_Ls)
    return calculate_appearance_fixation_time_from_clade_hmm(ts,fs,Ls)[0]

  
def calculate_fixed_weight(L,f):
    
    if L!=parse_file.EXTINCT:
        fixed_weight = f
    else:
        fixed_weight = 0
    
            
    return fixed_weight
    
def calculate_clade_fixed_weight(L,f):

    if L in parse_file.clade_fixed_states:
        fixed_weight = 1.0
    elif L in parse_file.clade_polymorphic_states:
        fixed_weight = 0.0
    else:
        fixed_weight = 0.0
        
    return fixed_weight
    