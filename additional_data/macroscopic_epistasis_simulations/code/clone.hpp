#ifndef CLONE_HPP
#define CLONE_HPP

#include <vector>

#include "object_pool.hpp"

typedef int Label; // used for identifying unique mutations
Label null_label = -1;

 
class LabelGenerator{        
    public:
        LabelGenerator(): current_label(0) {};
        Label get_next_label() {
            if(current_label == null_label){
                ++current_label;
            }
            return current_label++;
        };
        Label current_label;
};

class Mutation{
    public:
        Label label;
        double fitness_effect;
        int t;
};

Mutation null_mutation{null_label,1};

inline bool operator<(Mutation const & m1, Mutation const & m2){ return m1.label < m2.label; }  
inline bool operator==(Mutation const & m1, Mutation const & m2){ return m1.label == m2.label; }
inline bool operator!=(Mutation const & m1, Mutation const & m2){ return !(m1==m2); }

typedef std::vector<Mutation> MutationList;

const double WMAX = 1e100;

inline int sample_lineage_size(Random & random, double expected_size){
    // old version
    return sample_poisson(random, expected_size);
    // new version 
    //if(expected_size < 10){
    //    return sample_small_poisson(random, expected_size);
    //}
    //else{
    //    int size = sample_normal(random)*std::sqrt(expected_size)+expected_size-0.5;
    //    if(size < 0)
    //        size = 0;
    //    return size;
    //}
}

class Clone{
    public:
       int size;
       double fitness;
       double epistasis_parameter;
       MutationList mutations;

       Clone(){
           size=1;
           fitness=1;
           epistasis_parameter=-1;
       }

       void add_mutation(Mutation & mutation){
           fitness *= mutation.fitness_effect;
           if(fitness > WMAX){
               fitness = WMAX;
           }    
           epistasis_parameter = -1;
           mutations.push_back(mutation);
       }


       void add_mutation(Mutation && mutation){
           fitness *= mutation.fitness_effect;
           if(fitness > WMAX){
               fitness = WMAX;
           }    
           epistasis_parameter = -1;
           mutations.push_back(mutation);
       }


       void form_mutant(Clone const & clone, Mutation & mutation){
           size = 1.0;
           fitness = clone.fitness;
           mutations = clone.mutations;
           add_mutation(mutation);
       }

       void form_mutant(Clone const & clone, Mutation && mutation){
           size = 1.0;
           fitness = clone.fitness;
           mutations = clone.mutations;
           add_mutation(mutation);
       }
       
       Mutation & get_earliest_mutation(){
           if(mutations.size() > 0){
               return mutations.front();
            }
            else{
               return null_mutation;
            }
       }
       
       void remove_earliest_mutation(){
           mutations.erase(mutations.begin());
       }
};

typedef SharedObjectPool<Clone> ClonePool;
typedef std::vector<decltype(ClonePool().allocate())> Population;

#endif
