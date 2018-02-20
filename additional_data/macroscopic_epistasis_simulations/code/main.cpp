#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "stats.hpp" // to help with random number generation
#include "clone.hpp"
#include "dfe.hpp"

//#define DEBUG

typedef std::pair<Mutation, double> SNP; // A SNP is defined as a (mutation, frequency) pair
typedef std::vector<SNP> SNPList;

class Params{
    public:
        double U;
        double s;
        double g;
};

class Measurement{
    public:
        int time;
        double avg_fitness;
        std::vector<MutationList> clone_snps;
        SNPList snps;
};

typedef std::vector<Measurement> Results; // For returning results from an individual replicate population

typedef std::pair<double, double> ParsedMutation;

class ParsedResult {
    public:
        std::vector<double> xs;
        std::vector<double> ps;
        std::vector<ParsedMutation> ss;
};

typedef std::vector<ParsedResult> ParsedResultList;

std::vector<double> calculate_avg_xs(ParsedResultList const & parsed_results);
ParsedResult parse_result(Results const & results);
void print_results(Params const & params, ParsedResultList const & parsed_results);

template <class DFE>
Results evolve_population(Random & random, DFE dfe, double N, std::vector<int> times, int n, bool sequence_population=false);

std::vector<MutationList> sample_clones(Random & random, Population const & population, MutationList const & fixed_mutations, int n);
SNPList call_snps(Population const & population, MutationList const & fixed_mutations, const double threshold_frequency = 0.01);

int fitness_mcmc(int argc, char* argv[]);
int timecourse(int argc, char * argv[]);

const double dilution = 100;
const double Nmax = 5e+08; // maximum population size per well
const double N = Nmax/dilution; // population size at bottleneck
const double scale_factor = std::log2(dilution); // number of generations between bottlenecks
    
const int num_lines = 6; // number of compete, non-mutator lines
const int clone_sample_size=1; // number of clones sampled per time point
const int num_muts=5; // number of mutations measured in Ara-1

// sample times in the experiment (in days);
const std::vector<int> observed_ts{0, 75, 150, 225, 300, 375, 450, 525, 600, 675, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800, 1950, 2100, 2250, 2400, 2550, 2700, 2850, 3000, 3300, 3600, 3900, 4200, 4500, 4800, 5100, 5400, 5700, 6000, 6300, 6600, 6900, 7200, 7500};

// sampled fitnesses
const std::vector<double> observed_xs{0, 0.08923565, 0.14150781, 0.17350982, 0.16192306, 0.20482747, 0.20348269, 0.21211371, 0.22235515, 0.22377097, 0.20673246, 0.222472, 0.22523769, 0.24750729, 0.25841047, 0.23763903, 0.23185392, 0.24146727, 0.2619611,  0.2565393, 0.26906499, 0.26073441, 0.26448047, 0.29645801, 0.27273506,  0.26896062,  0.28210401,  0.28963463,  0.28350962,  0.28837608, 0.29721361, 0.31621516, 0.30258121, 0.31699151, 0.32151924, 0.31598268, 0.29787139,  0.31969785, 0.29828907, 0.31487821,  0.32589085}; 
const double experimental_sigma = 0.014;

int main(int argc, char * argv[]){
    // Simulation type is set here
    
    // MCMC algorithm
    //return fitness_mcmc(argc, argv);
    
    // Simulate mutation timecourse a la Lang et al (Nature, 2013)
    return timecourse(argc, argv);
}

int fitness_mcmc(int argc, char * argv[]){

    // Make sure there are enough command line arguments
    if(argc < 7){
        std::cout << "usage: " << argv[0] << " U0 s0 g0 dl num_runs num_steps\n";
        return 1;
    }

    // Read in parameters 
    double U0 = atof(argv[1]);
    double s0 = atof(argv[2]);
    double g0 = atof(argv[3]);
    double dl = atof(argv[4]);
    const int num_runs = atoi(argv[5]);
    const int num_steps = atoi(argv[6]);
    
    auto times = observed_ts;
    const int num_timepoints = times.size();

    // Create and seed random number generator
    Random random = create_random(0,true);
    auto draw_param_index = std::uniform_int_distribution<>(0,2); 
  
    double old_log_perror = -10000;
    double new_log_perror;
      
    Params old_params{-1,-1,-1};
    Params new_params{U0,s0,g0};
  
    for(int current_num_steps=0; current_num_steps < num_steps; ++current_num_steps){
        #ifdef DEBUG
            std::cerr << new_params.U << " " << new_params.s << " " << new_params.g << std::endl;
        #endif 
        
        // Create DFE
        // Delta DFE
        // auto dfe = create_tracked_dfe( DeltaDFE(new_params.U*scale_factor, new_params.s*scale_factor) );
        // Exponential DFE
        //auto dfe = create_tracked_dfe( ExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor) );        
        // Truncated Exponential DFE
        //auto dfe = create_tracked_dfe( TruncatedExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor, 4) );

        // Finite sites wierd DFE
        //const double eps = 3e-04;
        //auto dfe = create_finite_sites_dfe( random, WeirdDFE(new_params.U*scale_factor,new_params.s*scale_factor,eps), (int) new_params.g );
        //auto dfe = create_finite_sites_dfe( ExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor), (int) new_params.g);
        
        // Diminishing returns exponential
        auto dfe = create_tracked_dfe( create_wiser_dfe(ExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor), new_params.g / scale_factor) );
       
        // Diminishing returns delta
        //auto dfe = create_tracked_dfe( create_wiser_dfe(DeltaDFE(new_params.U*scale_factor, new_params.s*scale_factor), new_params.g / scale_factor) );

        //auto dfe = create_tracked_dfe( create_wiser_dfe( TruncatedExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor, 4), new_params.g / scale_factor ) );

        // peak exponential
        //auto dfe = create_tracked_dfe( create_sergey_dfe(ExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor), new_params.g / scale_factor) );

        //auto dfe = create_tracked_dfe( create_uncorrelated_dfe( ExponentialDFE(new_params.U*scale_factor, new_params.s*scale_factor) ) );
        
        ParsedResultList parsed_results;
        // measure first six fitness and mutation trajectories
        for(int i=0;i<num_lines;++i){
            #ifdef DEBUG
                std::cerr << i+1 << std::endl;
            #endif 
            // Evolve a replicate population
            auto results = evolve_population(random, dfe, N, times, clone_sample_size);
            parsed_results.push_back( parse_result( results ) );
        }
        
        // calculate the logarithm of the error probability 
        new_log_perror = 0;
        auto avg_xs = calculate_avg_xs(parsed_results);
        for(int i=1;i<num_timepoints;++i){
            new_log_perror -= (avg_xs[i]-observed_xs[i])*(avg_xs[i]-observed_xs[i])/2/experimental_sigma/experimental_sigma;
        } 
        
        if( (-sample_exponential(random) < new_log_perror - old_log_perror) || old_params.U < 0){
            // ACCEPTED!  
            #ifdef DEBUG
                std::cerr << "Accepted!" << std::endl;
            #endif 
            // update params
            old_params = new_params;
            old_log_perror = new_log_perror;
            // fill in the remaining populations 
            for(int i=num_lines;i<num_runs;++i){    
                #ifdef DEBUG
                    std::cerr << i+1 << std::endl;
                #endif 
                // Evolve a replicate population
                auto results = evolve_population(random, dfe, N, times, clone_sample_size);
                parsed_results.push_back( parse_result( results ) );
            }
            // print results
            print_results(old_params, parsed_results);    
        } 
        
        // generate new params
        new_params = old_params;
        switch(draw_param_index(random)){
            case 0:
                new_params.U *= std::exp( dl*sample_normal(random) );
                break;
            case 1:
                new_params.s *= std::exp( dl*sample_normal(random) );
                break;
            case 2:
                new_params.g *= std::exp( dl*sample_normal(random) );
                break;
            default:
                break;
        }
    }
    return 0;
} 

void print_results(Params const & params, ParsedResultList const & parsed_results){
    // params : xs1 ; ps1 ; ss1 : ... : xsn ; psn ; ssn
    std::cout << params.U << " " << params.s << " " << params.g;
    for(auto & parsed_result : parsed_results){
        std::cout << " : ";
        for(auto & x : parsed_result.xs){
            std::cout << x << " ";
        }
        std::cout << " ; ";
        for(auto & p : parsed_result.ps){
            std::cout << p << " ";
        }
        std::cout << " ; ";
        for(auto & parsed_mutation : parsed_result.ss){
            std::cout << parsed_mutation.first << "," << parsed_mutation.second << " ";
        }  
    }
    std::cout << std::endl;
}

int timecourse(int argc, char * argv[]){

    // Make sure there are enough command line arguments
    if(argc < 6){
        std::cout << "usage: " << argv[0] << " num_runs U0 s0 g0 Un\n";
        return 1;
    }

    const int num_runs = atoi(argv[1]);
    double U = atof(argv[2]);
    double s = atof(argv[3]);
    double g = atof(argv[4]);
    double Un = atof(argv[5]);
    const bool sequence_population = true;

    std::vector<int> times;
    for(int i=0;i<=120;++i){
        times.push_back(75*i);
    }
    const int num_timepoints = times.size();
    
    // print stuff to file
    std::cout << s << " " << U << " " << g << " " << Un << std::endl;
    // times
    for(auto & t : times){
        std::cout << (long) (t*scale_factor) << " ";
    }
    std::cout << std::endl;
    
    std::vector<double> xs(num_timepoints,0);
    std::vector<double> ps(num_timepoints,0);
    std::vector<double> ss(num_timepoints,0);

    // Create and seed random number generator
    Random random = create_random();
 
    // Create DFE
    
    //auto dfe = create_tracked_dfe( DeltaDFE(Ub,s0) );
    //auto dfe = create_tracked_dfe( create_composite_dfe( DeltaDFE(U*scale_factor,s*scale_factor), NeutralDFE(Un*scale_factor) ) );
    
    //auto dfe = create_tracked_dfe( create_composite_dfe( TruncatedExponentialDFE(U*scale_factor, s*scale_factor, 4), NeutralDFE(Un*scale_factor) ) );

    
    // Diminishing returns exponential DFE (Wiser et al, 2013)
    auto dfe = create_tracked_dfe( create_composite_dfe(create_wiser_dfe(ExponentialDFE(U*scale_factor, s*scale_factor), g / scale_factor), NeutralDFE(Un*scale_factor)) );

    for(int i=0;i<num_runs;++i){
        // Evolve a replicate population
        auto results = evolve_population(random, dfe, N, times, clone_sample_size, sequence_population);
        
        for(int j=0;j<num_timepoints;++j){
            auto & measurement = results[j];
            xs[j] = log(measurement.avg_fitness)/scale_factor;
            ps[j] = 0;
            for(auto & snps : measurement.clone_snps){
                /*if(snps.size() < 1){
                    std::cout << "No SNPS!\n";
                }
                else{
                    std::cout << snps.size() << std::endl;
                }*/
                ps[j] += snps.size();
            }
            ps[j] /= measurement.clone_snps.size();
        }
        
        // fitnesses
        for(auto & x : xs){
            std::cout << x << " ";
        }
        std::cout << ": ";
        
        // mutations
        for(auto & p : ps){
            std::cout << p << " ";
        }
        std::cout << ": ";
                
        // fitness effects along line of descent
        int j = 0;
        for(auto & mutation : results.back().clone_snps.front()){        
            std::cout << log(mutation.fitness_effect)/scale_factor << "," << mutation.t*scale_factor << " ";
            j += 1;
            if(j>5)
                break;
        }
        std::cout << ": ";
        std::cout << std::flush;
        // population metagenomic sequencing
        j=0;
        for(auto & measurement : results){
            j+=1;
            //std::cerr << j << " " << measurement.time << " " << measurement.snps.size() << std::endl;
            int k=0;
            for(auto & snp : measurement.snps){
                k+=1;
                //std::cerr << j << " " << k << std::endl;
                std::cout << snp.first.label << "," << snp.second << ";"; 
            }
            std::cout << " + ";
        }
        std::cout << "Done!";
        std::cout << std::endl;
    }
    return 0;
}


template <class DFE>
Results evolve_population(Random & random, DFE dfe, const double N, const std::vector<int> times, const int n, const bool sequence_population){

    Results results;
    MutationList fixed_mutations;
        
    ClonePool clone_pool(4*N*dfe.U0); // we use an object pool to make allocation of clones more efficient
    
    Population population; // the population in the current generation
    Population new_population; // the population in the next generation

    // initialize the population
    double total_fitness = 0;
    double total_size = 0;
    auto new_clone_ptr = clone_pool.allocate();
    new_clone_ptr->size = N;
    new_clone_ptr->fitness = 1.0;
    new_clone_ptr->mutations.clear();
    population.push_back(new_clone_ptr);
    total_size += new_clone_ptr->size;
    total_fitness += new_clone_ptr->size*new_clone_ptr->fitness;    

    auto observation_time_ptr = times.begin();
    auto observation_time_end = times.end();

    for(int t=0; observation_time_ptr < observation_time_end; ++t){
        
        
        if(t == *observation_time_ptr){
            // record stuff if is the right time
            
            std::cerr << t << std::endl;
        
            
            if(sequence_population){
                results.push_back(Measurement{t, total_fitness/total_size, sample_clones(random, population, fixed_mutations, n), call_snps(population, fixed_mutations)});
            }
            else{
                results.push_back(Measurement{t, total_fitness/total_size, sample_clones(random, population, fixed_mutations, n), SNPList()});
            }
            ++observation_time_ptr;
        }
        
        //std::cout << t << " " << total_size << " " << total_fitness << std::endl;
        double weighting_factor = N/total_fitness; // used to keep the total population size near N
        total_fitness = 0; // reset these variables so that we can 
        total_size = 0;    // remeasure them for the next generation

        for(auto & clone_ptr : population){

            //std::cout << "Drawing offspring...\n";
            // draw clonal and mutant offspring
            double U = dfe.get_mutation_rate(*clone_ptr);
            double expected_num_offspring = (clone_ptr->size)*(clone_ptr->fitness)*weighting_factor;
            int num_clones = sample_lineage_size(random, (1-U)*expected_num_offspring);
            int num_mutants = sample_lineage_size(random, U*expected_num_offspring);
            
            
            //std::cout << "Propagating clones...\n";
            if(num_clones > 0){
                // progagate clones
                clone_ptr->size = num_clones;
                new_population.push_back(clone_ptr);
                total_size += num_clones;
                total_fitness += num_clones*clone_ptr->fitness;
            }
            if(num_mutants > 0){
                // propagate mutants
                for(int i=0;i<num_mutants;++i){
                    new_clone_ptr = clone_pool.allocate();
                    new_clone_ptr->form_mutant(*clone_ptr, dfe.get_mutation(random, *clone_ptr,t));
                    new_population.push_back(new_clone_ptr);
                    total_fitness += new_clone_ptr->fitness;
                }
                total_size+=num_mutants;
            }
        }
        // make the next generation the current generation
        population.swap(new_population);
        new_population.clear();
                 
        // remove fixed mutations
        auto mutation = population.front()->get_earliest_mutation();
        if(mutation != null_mutation){
            bool fixed = true;
            for(auto & clone_ptr : population){
                if(clone_ptr->get_earliest_mutation() != mutation){
                    fixed = false;
                    break;
                }
            }
            if(fixed){
                for(auto & clone_ptr : population){
                    clone_ptr->remove_earliest_mutation();
                }
                dfe.fix_mutation(mutation);
                fixed_mutations.push_back(mutation);
            }
        }
                 
    }
    //std::cout << "Done! " << total_size << std::endl;
    return results;
}

std::vector<MutationList> sample_clones(Random & random, Population const & population, MutationList const & fixed_mutations, int n){

    std::vector<MutationList> clone_snps;
    double population_size=0;    
    for(auto const & clone_ptr : population){
            population_size += clone_ptr->size;
    }
    
    for(int i=0;i<n;++i){
        double location = sample_uniform(random)*population_size;
        double current_location = 0;
        for(auto const & clone_ptr : population){
            current_location += clone_ptr->size;
            if(location <= current_location){
                // found our clone!
                MutationList snps;
                // add the fixed mutations
                snps.insert(snps.end(),fixed_mutations.begin(),fixed_mutations.end());
                // add the additional clone mutations
                snps.insert(snps.end(), clone_ptr->mutations.begin(), clone_ptr->mutations.end());
                /*if(snps.size() < 1){
                    std::cout << "No SNPS!\n";
                }*/
                clone_snps.push_back(snps);
                break;
            }
        }
    }
    /*if(clone_snps.size() < n){
        std::cout << "Not enough clones!\n";
    }*/               
    return clone_snps;
}


SNPList call_snps(Population const & population, MutationList const & fixed_mutations, const double threshold_frequency){
    double N = 0;
    std::unordered_map<int,SNP> snp_map;

    for(auto const & clone_ptr : population){
        double clone_size = clone_ptr->size;
        N += clone_size;
        for(auto const & mutation : clone_ptr->mutations){
            if(!snp_map.count(mutation.label)){
                snp_map[mutation.label] = SNP{mutation, 0.0};
            }
            snp_map[mutation.label].second+=clone_size;
        }
    } 
    for(auto const & mutation : fixed_mutations){
        snp_map[mutation.label] = SNP{mutation, N};
    }
    
    // construct SNP List
    SNPList snps;
    snps.reserve(snp_map.size());
    for(auto & record : snp_map){
        record.second.second /= N;
        if(record.second.second >= threshold_frequency){
            snps.push_back(record.second);
        }
    }
    return snps;
}

std::vector<double> calculate_avg_xs(ParsedResultList const & parsed_results){
    std::vector<double> avg_xs; 
    for(int i=0, num_timepoints=parsed_results.front().xs.size(); i<num_timepoints;++i){
        avg_xs.push_back(0);
        for(auto & parsed_result : parsed_results){
            avg_xs[i] += parsed_result.xs[i];
        }
        avg_xs[i] /= parsed_results.size();
    }
    return avg_xs; 
}

ParsedResult parse_result(Results const & results){
    ParsedResult parsed_result;
    
    // fitnesses and mutation counts
    for(auto & measurement : results){
        double x = log(measurement.avg_fitness)/scale_factor;
        double p = 0;
        for(auto & snps : measurement.clone_snps){
            p+=snps.size();
        }
        p /= measurement.clone_snps.size();
        parsed_result.xs.push_back(x);
        parsed_result.ps.push_back(p);      
    }
    
    // fitness effects 
    int current_num_muts = 0;
    for(auto & mutation : results.back().clone_snps.front()){ 
        double s = log(mutation.fitness_effect)/scale_factor;
        double t = mutation.t*scale_factor;
        parsed_result.ss.push_back(ParsedMutation{s,t});
        current_num_muts += 1;
        if(current_num_muts>=num_muts)
                break;
    }
    
    return parsed_result;
}




