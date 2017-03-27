#include<vector>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <functional>
#include <algorithm>
#include <tuple>

#include "stats.hpp"
#include "trajectory.hpp"
#include "pvalue.hpp" 

std::tuple<double,double> calculate_combined_pvalue(Random & random, Trajectory const & trajectory, int max_num_bootstraps=10000, int min_num_bootstraps=10000);

std::tuple<int,double,double> calculate_deletion_pvalue(Random & random, Trajectory const & trajectory, int num_bootstraps);
std::tuple<int,double,double> calculate_duplication_pvalue(Random & random, Trajectory const & trajectory, int num_bootstraps);

const char delim=',';

int main(int argc, char * argv[]){

    // random number generator
    // deterministic seed to ensure reproducibility
    // once pipeline is completed
    auto random = create_random(42); 
    
    // used for reporting purposes
    int num_processed = 0;
    int num_passed = 0;
    int num_surprising = 0;
    
    // avg depth across the genome
    // used to trim trajectories with 
    // apparent deletions
    std::vector<double> avg_depths; 
    
    // iterate over all trajectory records in file
    std::string line;
    while(std::getline(std::cin,line)){
        
        // parse trajectory record
        std::stringstream line_stream(line);
        std::string item;
        std::string subitem;
        
        std::string allele;
        std::vector<double> times;
        std::vector<double> alts;
        std::vector<double> depths;
            
        // these entries not needed for this step in pipeline
        std::getline(line_stream, item, ','); // chromosome
        std::getline(line_stream, item, ','); // location
        std::getline(line_stream, allele, ','); // allele
            
        // times
        std::getline(line_stream, item, ',');
        std::stringstream item_stream0(item);
        while(std::getline(item_stream0, subitem, ' ')){
            if(subitem.size() > 0){
                times.push_back(std::stof(subitem));
            }
        }
        
        // alts
        std::getline(line_stream, item, ',');
        std::stringstream item_stream(item);
        while(std::getline(item_stream, subitem, ' ')){
            if(subitem.size() > 0){
                alts.push_back(std::stof(subitem));
            }
        }
            
        // depths
        std::getline(line_stream, item, ',');
        std::stringstream item_stream2(item);
        while(std::getline(item_stream2, subitem, ' ')){
            if(subitem.size() > 0){
                depths.push_back(std::stof(subitem));
            }
        }
        
        // the six statistics we will calculate
        // the first is whether it passed our basic filter
        bool passed_filter; 
        // the first three detect whether there is a deletion (and we should trim)
        int deletion_idx; 
        double fold_reduction, deletion_pvalue;
        // the next three detect whether there is a duplication (deal with this downstream)
        int duplication_idx; 
        double fold_increase, duplication_pvalue;
        // the next two are about the quality of the remaining trajectory
        double autocorrelation, combined_pvalue;       
        
        
        if(allele[1]=='D'){
            // Special trajectory containing average depths
            
            for(int i=0,imax=depths.size();i<imax;++i){
                if(times[i] > 100){
                    break; // don't include clone timespoints in here
                }
                avg_depths.push_back(depths[i]);
            }
            
            passed_filter = true;
            deletion_idx = 0;
            fold_reduction = 0;
            deletion_pvalue = 1;
            duplication_idx = 0;
            fold_increase = 0;
            duplication_pvalue = 1;
            autocorrelation = 0;
            combined_pvalue = 1;
            
        }
        else{
        
            // Regular trajectory
            
            // create trajectory
            auto trajectory = Trajectory();
            for(int i=0,imax=alts.size();i<imax;++i){
                if(times[i] > 100){
                    break; // don't include clone timepoints here
                }
                trajectory.push_back( Timepoint{ alts[i], depths[i], avg_depths[i] } );
            }
            
            num_processed+=1;
            
            if(num_processed % 1000 == 0){
                    std::cerr << num_processed << " trajectories processed, " << num_passed << " passed, " << num_surprising << " surprising!\n";
            }
            
            passed_filter = passes_filter(trajectory);
            
            if(!passed_filter){
                deletion_idx = 0;
                fold_reduction = 0;
                deletion_pvalue = 1;
                duplication_idx = 0;
                fold_increase = 0;
                duplication_pvalue = 1;
                autocorrelation = 0;
                combined_pvalue = 1;
            }
            else{
                
                if(avg_depths.size()==0){
                    // if this happens, trimming routine will choke
                    // should never get here?
                    deletion_idx = 0;
                    fold_reduction = 0;
                    deletion_pvalue = 1;
                    duplication_idx = 0;
                    fold_increase = 0;
                    duplication_pvalue = 1;
                }
                else{
                    
                    // test whether there is a deletion or duplication during the trajectory
                    std::tie(deletion_idx, fold_reduction, deletion_pvalue) = calculate_deletion_pvalue(random, trajectory, 10000);
                    std::tie(duplication_idx, fold_increase, duplication_pvalue) = calculate_duplication_pvalue(random, trajectory, 10000);
                    
                    // trim the trajectory at the deletion if necessary
                    if(deletion_pvalue < deletion_pvalue_threshold){
                        std::cerr << "Trimming trajectory: " << deletion_idx << " " << fold_reduction << " " << deletion_pvalue << std::endl;
                        for(int i=0,imax=trajectory.size()-deletion_idx;i<imax;++i){
                            trajectory.pop_back();
                        }
                        //std::cerr << "Done!" << std::endl;
                    }
                } // done with trimming
                
                
                // recalculate whether trajectory passes or not
                passed_filter = passes_filter(trajectory);
                
                
                if(!passed_filter){
                    autocorrelation = 0;
                    combined_pvalue = 1;
                }
                else{
                    // actually calculate autocorrelation score
                    num_passed+=1;
                    std::tie(autocorrelation, combined_pvalue) = calculate_combined_pvalue(random, trajectory); 
                }
            
                // increment number of "surprising" trajectories
                // used only for display purposes during run
                if(combined_pvalue <= 5e-02){
                    num_surprising += 1; 
                }
            } // done with pvalue calculations
        } // done with trajectory processing
        
        // print passed (but not necessarily significant) trajectory records 
        // (significance filtering happens later)
        if(passed_filter){
            std::cout << line << ", ";
            std::cout << autocorrelation << " " << combined_pvalue << ", ";
            std::cout << deletion_idx << " " << fold_reduction << " " << deletion_pvalue << ", ";
            std::cout << duplication_idx << " " << fold_increase << " " << duplication_pvalue << std::endl;  
        }  
    }
    std::cerr << "Finished: " << num_processed << " trajectories processed, " << num_passed << " passed, " << num_surprising << " surprising!\n";
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Calculates pvalue for trajectory
// based on several joint test statistics
//
// Returns: score, pvalue
//
////////////////////////////////////////////////////////////////////////////////


std::tuple<double,double> calculate_combined_pvalue(Random & random, Trajectory const & observed_trajectory, int max_num_bootstraps, int min_num_bootstraps){
     
    if(!passes_filter(observed_trajectory)){
        std::cerr << "Trajectory did not pass filter, should not get here!" << std::endl;
        std::tuple<double,double>{0,1};
    }
    
    // null model is based on permutation of observed trajectory
    // (with some added bells and whistles, see Supplementary Information)
    auto trajectory_generator = ResampledPermutationTrajectoryGenerator(observed_trajectory);
    
    // calculate num nonzero points
    int num_nonzero_points = 0;
    for(auto & timepoint : trajectory_generator.unmasked_observed_trajectory){
        if(timepoint.alt > 0){
            ++num_nonzero_points;
        }     
    }
    if(num_nonzero_points < 2){
        // can't have autocorrelation if only one timepoint with alt reads
        return std::tuple<double,double>{0, 1.0};   
    }

    class CombinedTestStatistic{
        public:
            double autocorrelation;
            double autocorrelation_pvalue;
            double max_run;
            double max_run_pvalue;
            double relaxation_time;
            double relaxation_time_pvalue;
            double combined_pvalue;
    };
    
    double observed_autocorrelation = calculate_autocorrelation(trajectory_generator.unmasked_observed_trajectory);
    double observed_max_run = calculate_max_run(trajectory_generator.unmasked_observed_trajectory);
    double observed_relaxation_time = calculate_relaxation_time(trajectory_generator.unmasked_observed_trajectory);
  
    CombinedTestStatistic observed_test_statistics = CombinedTestStatistic{observed_autocorrelation, 1.0, observed_max_run, 1.0, observed_relaxation_time, 1.0};
    
    std::vector<double> bootstrapped_autocorrelations;
    std::vector<double> bootstrapped_max_runs;
    std::vector<double> bootstrapped_relaxation_times;
    
    std::vector<CombinedTestStatistic> bootstrapped_test_statistics;
    
    //std::cout << observed_autocorrelation << " " << observed_max_run << " " << observed_relaxation_time << std::endl;
    
    // calculate # bootstrapped Ts > T
    int num_greater_autocorrelations = 0;
    int num_greater_max_runs = 0;
    int num_greater_relaxation_times = 0;
    int current_num_bootstraps = 0;    
    //std::cout << "Starting bootstraps!" << std::endl;
    while(current_num_bootstraps < max_num_bootstraps){
        
        current_num_bootstraps+=1;
        
        //if(current_num_bootstraps%10000==0){
            //std::cerr << current_num_bootstraps << std::endl;
        //}
        
        trajectory_generator.generate_bootstrapped_trajectory(random);
        
        // calculate bootstrapped test statistics
        double bootstrapped_autocorrelation = calculate_autocorrelation(trajectory_generator.unmasked_bootstrapped_trajectory);
        double bootstrapped_max_run = calculate_max_run(trajectory_generator.unmasked_bootstrapped_trajectory);      
        double bootstrapped_relaxation_time = calculate_relaxation_time(trajectory_generator.unmasked_bootstrapped_trajectory,true);
        
        // add them to corresponding vectors
        bootstrapped_autocorrelations.push_back(bootstrapped_autocorrelation);
        bootstrapped_max_runs.push_back(bootstrapped_max_run);
        bootstrapped_relaxation_times.push_back(bootstrapped_relaxation_time);

        // compare to observed values
        if( (bootstrapped_autocorrelation > observed_autocorrelation-1e-09) ){
            // bootstrapped trajectory is at least as extreme
            num_greater_autocorrelations+=1;
        }
        
        if( (bootstrapped_max_run > observed_max_run-1e-09) ){
            // bootstrapped trajectory is at least as extreme
            num_greater_max_runs+=1;
        }
        
        if( (bootstrapped_relaxation_time > observed_relaxation_time-1e-09) ){
            // bootstrapped trajectory is at least as extreme
            num_greater_relaxation_times+=1;
        }
        
        if( (num_greater_autocorrelations > min_numerator_counts) && (num_greater_max_runs > min_numerator_counts) && (num_greater_relaxation_times > min_numerator_counts) )
            break;
    }   
    
    // calculate pvalues
    double observed_autocorrelation_pvalue = calculate_pvalue_from_counts(num_greater_autocorrelations, current_num_bootstraps);
    double observed_max_run_pvalue = calculate_pvalue_from_counts(num_greater_max_runs, current_num_bootstraps);
    double observed_relaxation_time_pvalue = calculate_pvalue_from_counts(num_greater_relaxation_times, current_num_bootstraps);
    double observed_combined_pvalue = calculate_combined_pvalue(observed_autocorrelation_pvalue, observed_max_run_pvalue, observed_relaxation_time_pvalue);
    
    // short circuit if we know the clipped pvalue is already high
    if( observed_combined_pvalue == 1.0 ){
        return std::tuple<double, double>{observed_autocorrelation, 1.0};
    }
    
    //std::cerr << "Passed initial screen: " << observed_autocorrelation_pvalue << " " << observed_max_run_pvalue << " " << observed_relaxation_time_pvalue << " " << observed_combined_pvalue << std::endl;
    //std::cerr << "Generating auxilliary trajectories..." << std::endl;
    
    // otherwise, prepare to calculate bootstrapped combine pvalues
    // first populate at least ~ min num bootstrapped trajectories
    for(int i=0;i<min_num_bootstraps;++i){
    
        trajectory_generator.generate_bootstrapped_trajectory(random);
        
        // calculate bootstrapped test statistics
        double bootstrapped_autocorrelation = calculate_autocorrelation(trajectory_generator.unmasked_bootstrapped_trajectory);
        double bootstrapped_max_run = calculate_max_run(trajectory_generator.unmasked_bootstrapped_trajectory);      
        double bootstrapped_relaxation_time = calculate_relaxation_time(trajectory_generator.unmasked_bootstrapped_trajectory,true);
        
        // add them to corresponding vectors
        bootstrapped_autocorrelations.push_back(bootstrapped_autocorrelation);
        bootstrapped_max_runs.push_back(bootstrapped_max_run);
        bootstrapped_relaxation_times.push_back(bootstrapped_relaxation_time);
    }
    
    // then sort test statistic lists
    std::sort(bootstrapped_autocorrelations.begin(), bootstrapped_autocorrelations.end());
    std::sort(bootstrapped_max_runs.begin(), bootstrapped_max_runs.end());
    std::sort(bootstrapped_relaxation_times.begin(), bootstrapped_relaxation_times.end());
    
    // use these sorted lists to re-adjust observed pvalues
    //
    observed_autocorrelation_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_autocorrelations, observed_autocorrelation);
    //
    observed_max_run_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_max_runs, observed_max_run);
    //
    observed_relaxation_time_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_relaxation_times, observed_relaxation_time);
    //    
    observed_combined_pvalue = calculate_combined_pvalue(observed_autocorrelation_pvalue, observed_max_run_pvalue, observed_relaxation_time_pvalue);
    //
    
    //std::cerr << "Recalculated versions: " << observed_autocorrelation_pvalue << " " << observed_max_run_pvalue << " " << observed_relaxation_time_pvalue << " " << observed_combined_pvalue << std::endl;
    
    // now calculate a pvalue using the "combined pvalue" as a test statistic
    int num_greater = 0;
    current_num_bootstraps = 0;    
    while(current_num_bootstraps < max_num_bootstraps){
        
        current_num_bootstraps+=1;
        
        trajectory_generator.generate_bootstrapped_trajectory(random);
        
        // calculate bootstrapped test statistics
        double bootstrapped_autocorrelation = calculate_autocorrelation(trajectory_generator.unmasked_bootstrapped_trajectory);
        //
        double bootstrapped_max_run = calculate_max_run(trajectory_generator.unmasked_bootstrapped_trajectory); 
        //     
        double bootstrapped_relaxation_time = calculate_relaxation_time(trajectory_generator.unmasked_bootstrapped_trajectory,true);
        //
        
        // calculate pvalues for those test statistics
        double bootstrapped_autocorrelation_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_autocorrelations, bootstrapped_autocorrelation);
        //
        double bootstrapped_max_run_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_max_runs, bootstrapped_max_run);
        //
        double bootstrapped_relaxation_time_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_relaxation_times, bootstrapped_relaxation_time);
        //
        double bootstrapped_combined_pvalue = calculate_combined_pvalue(bootstrapped_autocorrelation_pvalue, bootstrapped_max_run_pvalue, bootstrapped_relaxation_time_pvalue);
        //
          
        //std::cerr << bootstrapped_autocorrelation_pvalue << " " << bootstrapped_max_run_pvalue << " " << bootstrapped_relaxation_time_pvalue << " " << bootstrapped_combined_pvalue << std::endl;
          
        if(bootstrapped_combined_pvalue <= observed_combined_pvalue){
            num_greater += 1;
        }    
    
        if(num_greater > min_numerator_counts){
            break;
        }
    }
    
    // calculate final pvalue
    double pvalue = calculate_pvalue_from_counts(num_greater, current_num_bootstraps);
    
    // decrease uncertainty for things near the boundary
    if(pvalue>0.01 && pvalue<0.05){
    
        // run it again
        // now calculate a pvalue using the "combined pvalue" as a test statistic
        while(current_num_bootstraps < 1000000){
        
            current_num_bootstraps+=1;
        
            trajectory_generator.generate_bootstrapped_trajectory(random);
        
            // calculate bootstrapped test statistics
            double bootstrapped_autocorrelation = calculate_autocorrelation(trajectory_generator.unmasked_bootstrapped_trajectory);
            //
            double bootstrapped_max_run = calculate_max_run(trajectory_generator.unmasked_bootstrapped_trajectory); 
            //     
            double bootstrapped_relaxation_time = calculate_relaxation_time(trajectory_generator.unmasked_bootstrapped_trajectory,true);
            //
        
            // calculate pvalues for those test statistics
            double bootstrapped_autocorrelation_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_autocorrelations, bootstrapped_autocorrelation);
            //
            double bootstrapped_max_run_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_max_runs, bootstrapped_max_run);
            //
            double bootstrapped_relaxation_time_pvalue = calculate_pvalue_from_sorted_list(bootstrapped_relaxation_times, bootstrapped_relaxation_time);
            //
            double bootstrapped_combined_pvalue = calculate_combined_pvalue(bootstrapped_autocorrelation_pvalue, bootstrapped_max_run_pvalue, bootstrapped_relaxation_time_pvalue);
            //
          
            if(bootstrapped_combined_pvalue <= observed_combined_pvalue){
                num_greater += 1;
            }    
    
        }
    
        // calculate final pvalue
        pvalue = calculate_pvalue_from_counts(num_greater, current_num_bootstraps);
    
    }
    
    //std::cerr << "Final pvalue: " << pvalue << std::endl;
    
    return std::tuple<double,double>{observed_autocorrelation, pvalue};
}

////////////////////////////////////////////////////////////////////////////////
//
// Calculates pvalue and location for trimming trajectory
// if it looks like a deletion has occured
//
// Returns: index, fold_change, pvalue
//
////////////////////////////////////////////////////////////////////////////////
std::tuple<int,double,double> calculate_deletion_pvalue(Random & random, Trajectory const & trajectory, int num_bootstraps){

    std::vector<double> depth_trajectory;
    std::vector<int> ks;

    // construct depth trajectory
    for(int k=0,kmax=trajectory.size();k<kmax;++k){
        if(trajectory[k].median_depth < min_depth){
            // do nothing!
        }
        else{
            depth_trajectory.push_back(trajectory[k].depth/trajectory[k].median_depth);
            ks.push_back(k);
        }
    }
    
    int L = depth_trajectory.size();
    
    if(L < 1){
        return std::tuple<int,double,double>{-1,0,1};
    }
    
    int observed_k;
    double observed_fold_change, observed_loglikelihood;
    int dummy_observed_k;
    double dummy_observed_fold_change;
    std::tie(observed_k, observed_fold_change, observed_loglikelihood) = calculate_depth_change(depth_trajectory);
    if(observed_fold_change > deletion_fold_reduction_threshold){
        std::tie(observed_k, observed_fold_change, observed_loglikelihood) = calculate_depth_change(depth_trajectory, deletion_fold_reduction_threshold);    
    }
    
    int cutoff_idx = ks[L-observed_k];
    
    if(cutoff_idx > trajectory.size()){
        std::cerr << "Panic! " << observed_k << " " << cutoff_idx << " " << trajectory.size() << std::endl;
    }
    int num_greater = 0;
    int current_num_bootstraps = 0;
    for(;current_num_bootstraps < num_bootstraps; ++current_num_bootstraps){
    
        if(num_greater > 10){
            break;
        }
    
        std::shuffle(depth_trajectory.begin(), depth_trajectory.end(), random);
        int bootstrapped_k;
        double bootstrapped_fold_change, bootstrapped_loglikelihood;
        int dummy_bootstrapped_k;
        double dummy_bootstrapped_fold_change;
        std::tie(bootstrapped_k, bootstrapped_fold_change, bootstrapped_loglikelihood) = calculate_depth_change(depth_trajectory);
        if(bootstrapped_fold_change > deletion_fold_reduction_threshold){
            std::tie(bootstrapped_k, dummy_bootstrapped_fold_change, bootstrapped_loglikelihood) = calculate_depth_change(depth_trajectory, deletion_fold_reduction_threshold); 
        }
        
        if(bootstrapped_loglikelihood >= observed_loglikelihood){
            num_greater += 1;
        }
    }
    
    double pvalue = (num_greater+1)*1.0/(current_num_bootstraps+1);
    return std::tuple<int,double,double>{cutoff_idx, observed_fold_change, pvalue};
}


////////////////////////////////////////////////////////////////////////////////
//
// Calculates pvalue and location for trimming trajectory
// if it looks like a duplication has occured
//
// Returns: index, fold_change, pvalue
//
////////////////////////////////////////////////////////////////////////////////
std::tuple<int,double,double> calculate_duplication_pvalue(Random & random, Trajectory const & trajectory, int num_bootstraps){

    std::vector<double> depth_trajectory;
    std::vector<int> ks;

    // construct depth trajectory
    for(int k=0,kmax=trajectory.size();k<kmax;++k){
        if(trajectory[k].median_depth < min_depth){
            // do nothing!
        }
        else{
            depth_trajectory.push_back(trajectory[k].depth/trajectory[k].median_depth);
            ks.push_back(k);
        }
    }
    
    int L = depth_trajectory.size();
    
    if(L < 1){
        return std::tuple<int,double,double>{-1,0,1};
    }
    
    int observed_k;
    double observed_fold_change, observed_loglikelihood;
    int dummy_observed_k;
    double dummy_observed_fold_change;
    std::tie(observed_k, observed_fold_change, observed_loglikelihood) = calculate_depth_change(depth_trajectory);
    if(observed_fold_change < duplication_fold_increase_threshold){
        std::tie(observed_k, observed_fold_change, observed_loglikelihood) = calculate_depth_change(depth_trajectory, duplication_fold_increase_threshold);    
    }
    
    int cutoff_idx = ks[L-observed_k];
    
    if(cutoff_idx > trajectory.size()){
        std::cerr << "Panic! " << observed_k << " " << cutoff_idx << " " << trajectory.size() << std::endl;
    }
    int num_greater = 0;
    int current_num_bootstraps = 0;
    for(;current_num_bootstraps < num_bootstraps; ++current_num_bootstraps){
    
        if(num_greater > 10){
            break;
        }
    
        std::shuffle(depth_trajectory.begin(), depth_trajectory.end(), random);
        int bootstrapped_k;
        double bootstrapped_fold_change, bootstrapped_loglikelihood;
        int dummy_bootstrapped_k;
        double dummy_bootstrapped_fold_change;
        std::tie(bootstrapped_k, bootstrapped_fold_change, bootstrapped_loglikelihood) = calculate_depth_change(depth_trajectory);
        if(bootstrapped_fold_change < duplication_fold_increase_threshold){
            std::tie(bootstrapped_k, bootstrapped_fold_change, bootstrapped_loglikelihood) = calculate_depth_change(depth_trajectory, duplication_fold_increase_threshold); 
        }
        
        if(bootstrapped_loglikelihood >= observed_loglikelihood){
            num_greater += 1;
        }
    }
    
    double pvalue = (num_greater+1)*1.0/(current_num_bootstraps+1);
    return std::tuple<int,double,double>{cutoff_idx, observed_fold_change, pvalue};
}
