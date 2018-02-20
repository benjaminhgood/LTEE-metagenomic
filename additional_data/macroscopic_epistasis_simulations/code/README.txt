LICENSE: The code located here is provided as-is, with no warranty, etc. (It's under GPL v2.) 

REQUIREMENTS: Python (v2.7) with the NumPy (v1.7) and SciPy (v0.12) libraries, and a C++ compiler that supports the C++11 standard (e.g., GNU's g++). 

USAGE: 

POPULATION GENETIC SIMULATIONS

The simulation type and DFE are set in the main.cpp file at the marked locations. The simulations can be compiled with the command "g++ -std=c++11 *.cpp -O3 -o simulate_ltee". The default setting is a Markov Chain Monte Carlo (MCMC) algorithm with the fitness-mediated epistasis model proposed by Wiser et al (Science, 2013). The simulation is run using the command

./simulate_ltee Ub sigma g 0.1 num_lines 1

where Ub is the beneficial mutation rate, sigma is the mean of the exponential DFE, g=1/(2*Xc) controls the strength of the diminishing returns epistasis, and num_lines is the number of independent lines to simulate. The output of the simulation is printed to stdout. It is of the form:

Ub, sb, g : X1_t1 X1_t2 ... X1_tmax ; M1_t1 M1_t2 ... M1_tmax; s1_1,t1_1 s1_2,t1_2 ... s1_5,t1_5 : X2_t1 X2_t2 ...

Here Xi_tk is the average fitness of replicate population i at timepoint k, Mi_tk is the number of mutations in a clone sampled from replicate population i at timepoint k, si_k is the fitness effect of the kth fixed mutation in replicate population i, and ti_k is the origination time of the kth fixed mutation in replicate population i. 

RUNTIME: Execution time strongly depends on the DFE parameters (and the mutation rate in particular), since these dictate how many independent lineages exist in the population. Typical runtimes on a laptop computer range from < 1s per replicate population for Ub ~ 1e-07 to 1hr per replicate for Ub ~ 1e-04. 

APPROXIMATE LIKELIHOOD CALCULATION

We have also included a reference implementation of our bootstrapped approximate likelihood calculation located in the file "fit.py". The function is called as

calculate_simulation_loglikelihood_pvalue(simulated_trajectories, observed_trajectory, sigma_measurement, simulated_dxs, observed_dx, sigma_dx_measurement, intercept=None, num_bootstraps = 10000)

The parameters are:

simulated_trajectories --- a (n,k) numpy array of fitness values whose k columns are timepoints and whose n rows are simulated replicate populations

observed_trajectory --- a (1,k) numpy array of fitness values whose k columns are timepoints for the observed fitness trajectory

sigma_measurement --- the estimated measurement error for the fitness measurements

simulated_dxs --- a (n,1) numpy array of fitness changes between generation 40,000 and 50,000 whose n rows are simulated replicate populations. 

observed_dx --- the observed change in fitness between generation 40,000 and 50,000 

sigma_dx_measurement --- the estimated measurement error for the change in fitness measurements. 

intercept --- the fitness at time 0. If None, the algorithm attempts to fit the intercept numerically using maximum likelihood. 

num_bootstraps --- number of bootstrap samples to simulate 

The output is a tuple (loglikelihood, pvalue, b), where the first item is the estimated loglikelihood, the second item is the pvalue for the chisquared goodness of fit test, and b is the estimated intercept.  

JOURNAL REFERENCE: Good BH and Desai MM, "The impact of macroscopic epistasis on long-term evolutionary dynamics," arXiv:1408.4108 (2014). 
