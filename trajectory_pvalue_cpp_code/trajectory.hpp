#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include<vector>
#include <iostream>
#include <random>
#include <functional>
#include <algorithm>
#include <tuple>

#include "stats.hpp"

const double min_depth = 5;
const double min_initial_depth = 10;

const double autocorrelation_threshold = 0.2;

const double deletion_pvalue_threshold = 1e-02;
const double deletion_fold_reduction_threshold = -1.0;
const double duplication_fold_increase_threshold = 0.9;

class Timepoint{
    public:
        double alt;
        double depth;
        double median_depth;
};

typedef std::vector<Timepoint> Trajectory;

inline bool filter_timepoint(Timepoint const & timepoint){
    if((timepoint.median_depth < min_depth) || (timepoint.depth < min_depth)){
        return true;
    }
    else{
        return false;
    } 
}

/* Filters out trajectories that do not have at least one timepoint with 
   2 or more alternate alleles */
inline bool passes_filter(Trajectory const & trajectory){

    if(trajectory[0].depth < min_initial_depth){
        return false;
    }

    double initial_frequency = trajectory[0].alt/trajectory[0].depth;
    
    if(initial_frequency>0.1){
        return false;
    }
    
    double max_frequency = initial_frequency;
    double min_frequency = initial_frequency;
    double avg_frequency = 0;

    int num_timepoints = 0; // # timepoints that pass depth filter
    int num_confirmed_timepoints = 0; // # of timepoints w/ >= 2 alts
    int num_good_timepoints = 0; // # timepoints w/ >= 2 alts and >= 5% freq
    int num_zero_timepoints = 0; // # timepoints w/ < 1 alts
    
    for(int i=1,imax=trajectory.size();i<imax;++i){
        auto const & timepoint = trajectory[i];
        if(!filter_timepoint(timepoint)){
        
            double frequency = timepoint.alt/timepoint.depth;

            // update average max and min            
            avg_frequency += frequency;
            num_timepoints += 1;

            max_frequency = (frequency > max_frequency) ? frequency : max_frequency;
            min_frequency = (frequency < min_frequency) ? frequency : min_frequency;
            
            if (timepoint.alt < 1){
                num_zero_timepoints += 1;
            }
            
            if(timepoint.alt >= 2){
                num_confirmed_timepoints += 1;
            }
            
            if(timepoint.alt >= 3 && frequency >= 0.10){
                num_good_timepoints += 1;
            }
        }
    }
    
    avg_frequency /= num_timepoints;
    // cap average at 50%
    avg_frequency = (avg_frequency > 0.5) ? 0.5 : avg_frequency;
    
    // ensure that there are at least some alts to work with
    
    if(num_good_timepoints < 1){
        return false;
    }
    
    if(num_confirmed_timepoints < 3){
        return false;
    }
    
    // weed out things that hang out at ~ 10% frequency forever
    //
    //     not at zero very much        max change isn't very high
    if( ((max_frequency-avg_frequency) < 0.10) ){
        return false;
    }
    
    // passes all the filters, so let's go! 
    return true;
    
}

class FrequencyGenerator{
    public:
        FrequencyGenerator(): sample_x(0), sample_y(0) {};
        FrequencyGenerator(double p, double D): p(p), D(D), sample_x(p*D), sample_y((1-p)*D){};
        
        double sample_p(Random & random){
            return p; 
            
            /*double x = sample_x(random);
            double y = sample_y(random);
            return x/(x+y);*/
        };
        
    double p;
    double D;
    
    private:
        std::gamma_distribution<double> sample_x, sample_y;
};

class ResampledPermutationTrajectoryGenerator{
    public:
    
        Trajectory bootstrapped_trajectory; // a full bootstrapped replica of the observed trajectory
        Trajectory unmasked_observed_trajectory; // only non-filtered points in the observed trajectory
        Trajectory unmasked_bootstrapped_trajectory; // only non-filtered points in the bootstrapped trajectory 
        
        double avg_frequency;
        
        ResampledPermutationTrajectoryGenerator(Trajectory const & trajectory){
        
            double total_alts = 0;
            double total_depths = 0;
            double frequency_multiplier = 1;
        
            // calculate avg frequency
            for(int i=0,imax=trajectory.size();i<imax;++i){
                Timepoint const & timepoint = trajectory[i]; 
                if(!filter_timepoint(timepoint)){
                    total_alts += timepoint.alt;
                    total_depths += timepoint.depth;
                }   
            }
            
            avg_frequency = total_alts/total_depths;
            
            if(avg_frequency > 0.5){
                frequency_multiplier = 0.5/avg_frequency;
            }
        
            for(int i=0,imax=trajectory.size();i<imax;++i){
                Timepoint const & timepoint = trajectory[i];
                bootstrapped_trajectory.push_back(timepoint);
                if(!filter_timepoint(timepoint)){
                    
                    unmasked_indices.push_back(i);
                    unmasked_observed_trajectory.push_back(timepoint);
                    unmasked_bootstrapped_trajectory.push_back(timepoint);
                    
                    double frequency = timepoint.alt/timepoint.depth*frequency_multiplier;
                    frequency = (frequency > 1e-03) ? frequency : 1e-03;
                    
                    ps.push_back(FrequencyGenerator(frequency,timepoint.depth));
                }
            }
            
            for(int i=0,imax=ps.size();i<imax;++i){
                permuted_p_indices.push_back(i);
            } 
            
            
        };
        
        Trajectory & generate_bootstrapped_trajectory(Random & random){
            do{
                // shuffle estimated ps
                std::shuffle(permuted_p_indices.begin()+1,permuted_p_indices.end(),random);
                
                // resample alts 
                for(int i=1,imax=unmasked_bootstrapped_trajectory.size();i<imax;++i){
                    //double A = round(ps[i]*unmasked_bootstrapped_trajectory[i].depth);
                    
                    double p = ps[permuted_p_indices[i]].sample_p(random);
                    double D = unmasked_bootstrapped_trajectory[i].depth;
                    double A = sample_poisson(random, D*p);
                    double R = sample_poisson(random, D*(1-p));
                    A = round(A/(A+R)*D);
                    unmasked_bootstrapped_trajectory[i].alt = A;
                    bootstrapped_trajectory[unmasked_indices[i]].alt = A;    
                }
            }while(!passes_filter(bootstrapped_trajectory));
            return bootstrapped_trajectory;
        };
      
    private:
        std::vector<int> unmasked_indices;
        std::vector<FrequencyGenerator> ps; 
        std::vector<int> permuted_p_indices;
};


/*inline double calculate_max_run(Trajectory const & trajectory){
    double max_run = 0;
    double current_run = 0;
    for(auto const & timepoint : trajectory){
        double frequency = timepoint.alt/timepoint.depth;
        if(frequency > 0){
            current_run += timepoint.alt/timepoint.depth;
        }
        else{
            max_run = (current_run > max_run) ? current_run : max_run;
            current_run = 0;
        }
    }
    return max_run;
}*/

inline double calculate_max_run(Trajectory const & trajectory){

    // calculate average frequency;
    double total_alts = 0;
    double total_depths = 0;
    double num_zero_timepoints = 0;
    double num_timepoints = 0;
    for(auto const & timepoint : trajectory){
        total_alts += timepoint.alt;
        total_depths += timepoint.depth;
        if(total_alts < 1){
            num_zero_timepoints += 1;
        }   
        num_timepoints += 1;
    }
    double avg_frequency = total_alts/total_depths;
    // cap it at a certain value
    avg_frequency = (avg_frequency < 0.5) ? avg_frequency : 0.5;
    
    if(num_zero_timepoints > 0.3*num_timepoints){
        // zero is close to median, weight avg toward that
        avg_frequency = avg_frequency*exp(-(num_zero_timepoints-0.3*num_timepoints)/5)/(1+exp(-(num_zero_timepoints-0.3*num_timepoints)/5));
    }
    else{
        // mostly nonzero timepoints, use avg
        avg_frequency = avg_frequency/(1+exp(-(0.3*num_timepoints-num_zero_timepoints)/5));
    }
    
    // if enough zero points... floor is zero. 
    // otherwise, floor is average. 
    
    //std::cout << "avg freq=" << avg_frequency << std::endl;
    
    double max_run = 0;
    double current_run = 0;
    double current_run_length = 0;
    for(auto const & timepoint : trajectory){
        double f = timepoint.alt/timepoint.depth;
        if(timepoint.alt - avg_frequency*timepoint.depth >= 1){
            // above by at least a read! 
            //double d_run = (f-avg_frequency)*timepoint.depth;
            //d_run = (d_run < 1) ? 0 : d_run;
            //current_run += d_run/timepoint.depth;
            current_run += (f-avg_frequency);
            current_run_length += 1;
            //std::cout << f << " " << avg_frequency << " " << current_run << std::endl;
        }
        else{
            //std::cout << "Ending! " << current_run << std::endl;
            max_run = (current_run > max_run && current_run_length>1) ? current_run : max_run;
            current_run = 0;
            current_run_length = 0;
        }
    }
    max_run = (current_run > max_run && current_run_length>1) ? current_run : max_run;
    return max_run;
}

inline double calculate_autocorrelation(Trajectory const & trajectory){
    
    // weighted autocorrelation, weighted average, etc. 
    
    // first calculate average
    double avg_f = 0;
    double average_weight = 0;
    for(auto & timepoint : trajectory){
        double f = timepoint.alt/timepoint.depth;
        double weight = timepoint.depth;
        avg_f += f*weight;
        average_weight += weight;
    }
    avg_f /= average_weight;
    
    
    // cap the avg frequency at 0.5
    double p = (avg_f > 0.5) ? 0.5 : avg_f; 
    
    double autocovariance=0;
    double covariance_weight=0;
    for(int i=1,imax=trajectory.size();i<imax;++i){
        double f = trajectory[i].alt/trajectory[i].depth;
        double fprev = trajectory[i-1].alt/trajectory[i-1].depth;
        double weight = sqrt(trajectory[i].depth*trajectory[i-1].depth);
        double read_difference = trajectory[i].alt-trajectory[i].depth*p;
        double read_difference_prev=  trajectory[i-1].alt-trajectory[i-1].depth*p;
        if ((read_difference >= 1 || read_difference <= -1) && (read_difference_prev >= 1 || read_difference_prev <= -1)){
            autocovariance += (f-p)*(fprev-p)*weight;
        }
        covariance_weight += weight;
    }
    autocovariance /= covariance_weight;
    
    
    double variance=0;
    double variance_weight=0;
    for(int i=0,imax=trajectory.size();i<imax;++i){
        double f = trajectory[i].alt/trajectory[i].depth;
        double weight = trajectory[i].depth;
        variance += (f-p)*(f-p)*weight;
        variance_weight += weight;
    }
    variance /= variance_weight;
    
    double autocorrelation = autocovariance/p/p;
    autocorrelation = (autocorrelation < 0) ? 0 : autocorrelation;
    
    //std::cout << "Avg f = " << avg_f << " C = " << autocorrelation << " sigma = " << sqrt(variance) << std::endl;
       
    
    return sqrt(autocorrelation);
    
}

inline double calculate_relaxation_time(Trajectory const & trajectory, bool is_bootstrap=false){
    
    // weighted autocorrelation, weighted average, etc. 
    
    // first calculate average
    double avg_f = 0;
    double average_weight = 0;
    double num_zero_timepoints = 0;
    double num_timepoints = 0;
    for(auto & timepoint : trajectory){
        double f = timepoint.alt/timepoint.depth;
        double weight = timepoint.depth;
        avg_f += f*weight;
        average_weight += weight;
        if(timepoint.alt < 1){
            num_zero_timepoints+=1;
        }
        num_timepoints += 1;
    }
    avg_f /= average_weight;
    
    //std::cout << avg_f << std::endl;
    
    if(num_zero_timepoints > 0.3*num_timepoints){
        if(is_bootstrap){
            return trajectory.size();
        }
        else{
            return -1;
        }
    }
    
    double running_avg_f = 0;
    double running_avg_weight = 0;
    double relaxation_time = -1;
    for(int i=1,imax=trajectory.size();i<imax;++i){
        double f = trajectory[i].alt/trajectory[i].depth;
        double weight = 1.0;
        running_avg_f += f*weight;
        running_avg_weight += weight;
        
        //std::cout << running_avg_f/running_avg_weight << std::endl;
        
        if( (i>5) and (running_avg_f/running_avg_weight > avg_f*0.60) ){
            relaxation_time = i;
            break;
        }
        
    }
    
    // should rarely happen
    if(relaxation_time < 0){
        
        relaxation_time = trajectory.size();
    }
    
    
    return relaxation_time;
      
}

inline std::tuple<int, double, double> calculate_depth_change(std::vector<double> & depth_trajectory, double fold_change=0){

    int L = depth_trajectory.size();
    double m1 = 0;
    double m21 = 0;
    double m2 = 0;
    double m22 = 0;

    // Calculate total variance
    for(auto const & x : depth_trajectory){
        m1 += x;
        m21 += x*x;
    }
    m1/=L;
    m21/=L;
    
    double var0 = m21 - m1*m1;
    
    double max_loglikelihood=0;
    int max_k=0;
    double max_fold_change=0;
    bool first=true;
    for(int k=1;k<L;++k){
        double x = depth_trajectory[L-k];
        m1 = ((L-k+1)*m1 - x)/(L-k);
        m21 = ((L-k+1)*m21 - x * x)/(L-k);
        m2 = ((k-1)*m2 + x)/k;
        m22 = ((k-1)*m22 + x * x)/k;
        
        double var1 = (m21-m1*m1);
        double mavg;
        if(fold_change==0){
            mavg=m2;
        }
        else{
            mavg=m1*std::pow(2,fold_change);
        }
        double var2 = (m22-2*mavg*m2+m2*m2);
        
        double loglikelihood = L*std::log( L * var0 / ( (L-k) * var1 + k * var2 ) ); 
        
        if(loglikelihood > max_loglikelihood || first){
            first=false;
            max_loglikelihood = loglikelihood;
            max_k = k;
            if(m1 > 0 && m2 > 0){
                max_fold_change = std::log2(mavg/m1);
            }
            else{
                max_fold_change = 1e10;
            }
        }
    }
    
    return std::tuple<int,double,double>{max_k, max_fold_change, max_loglikelihood};

}

#endif