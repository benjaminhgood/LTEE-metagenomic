#ifndef DFE_HPP
#define DFE_HPP

#include "stats.hpp"
#include "clone.hpp"
#include <cmath>
#include <memory>

const double UMAX = 0.01;

// A dummy DFE that only returns neutral (i.e., s=0 or W=1) fitness effects
class NeutralDFE{
    public:
        const double U0;
        NeutralDFE(double U): U0(U) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const { return 1; };
};

// A DFE where all mutations share the same fitness effect 
class DeltaDFE{
    public:
        const double U0;
        const double W;
        DeltaDFE(double U, double s): U0(U), W(std::exp(s)) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const { return W; } ;
};

// An exponentially distributed DFE with an average fitness effect savg
class ExponentialDFE{
    public:
        const double U0;
        const double W;
        ExponentialDFE(double U, double savg): U0(U), W(std::exp(savg)) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const { 
            return std::pow(W,sample_exponential(random)); 
        };
};

// A uniformly distributed DFE on the interval [0,smax]
class UniformDFE{
    public:
        const double U0;
        const double W;
        UniformDFE(double U, double smax): U0(U), W(std::exp(smax)) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return std::pow(W,sample_uniform(random));
        };
};

// A truncated exponential DFE 
class TruncatedExponentialDFE{
    public:
        const double U0;
        const double W;
        const double p_thresh;
        TruncatedExponentialDFE(double U, double smax, double c): U0(U), W(std::exp(smax)), p_thresh(1-std::exp(-c)) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return std::pow(W,-log(1.0-sample_uniform(random)*p_thresh));
        };
};
      

// A DFE predicted to have Wiser-like fitness trajectories in the weak mutation limit
class WeirdDFE{
    public:
        const double U0;
        const double W;
        const double eps;
        WeirdDFE(double U, double sigma, double eps): U0(U), W(std::exp(sigma)), eps(eps) {};
        double get_mutation_rate(Clone const & individual) const { return U0; };
        double get_fitness_effect(Random & random, Clone const & individual) const {
            double x;
            
            do{ x = eps/sample_uniform(random); }  // draw pareto variable
            while( sample_uniform(random) > std::exp(x-eps) ); // rejection sample
            
            return std::pow(W,x);
        };
};

template <class DFE1, class DFE2>
class CompositeDFE{
    public:
        DFE1 dfe1;
        DFE2 dfe2;
        double U0;
        CompositeDFE(DFE1 dfe1, DFE2 dfe2): dfe1(dfe1), dfe2(dfe2) {};
        double get_mutation_rate(Clone const & individual) const { return dfe1.get_mutation_rate(individual)+dfe2.get_mutation_rate(individual); };
        double get_fitness_effect(Random & random, Clone const & individual){
            if( sample_uniform(random) < dfe1.get_mutation_rate(individual)/(dfe1.get_mutation_rate(individual)+dfe2.get_mutation_rate(individual)) ){
                return dfe1.get_fitness_effect(random,individual);
            }
            else{
                return dfe2.get_fitness_effect(random,individual);
            }
        }
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE1, class DFE2>
inline auto create_composite_dfe(DFE1 dfe1, DFE2 dfe2) -> CompositeDFE<DFE1, DFE2> {
    return CompositeDFE<DFE1, DFE2>(dfe1, dfe2);
}

// Turn a DFE into an uncorrelated DFE, i.e., rho(s) = f(X+s)
template <class DFE>
class UncorrelatedDFE{
    public:
        DFE dfe;
        const double U0;
        UncorrelatedDFE(DFE dfe): dfe(dfe), U0(dfe.U0) {};
        double get_mutation_rate(Clone const & individual) const { return dfe.get_mutation_rate(individual); };
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return dfe.get_fitness_effect(random, individual)/individual.fitness;
        }; 
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_uncorrelated_dfe(DFE dfe) -> UncorrelatedDFE<DFE> {
    return UncorrelatedDFE<DFE>(dfe);
}

// Turn a DFE into a Wiser DFE, i.e., s(X) = s(0)*exp(-g*X)
template <class DFE>
class WiserDFE{
    public:
        DFE dfe;
        const double U0;
        const double g; 
        WiserDFE(DFE dfe, double g): dfe(dfe), U0(dfe.U0), g(g) {};       
        double get_mutation_rate(Clone const & individual) const { return dfe.get_mutation_rate(individual); };
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return std::pow(dfe.get_fitness_effect(random, individual), std::pow(individual.fitness,-g));
        };
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_wiser_dfe(DFE dfe, double g) -> WiserDFE<DFE> {
    return WiserDFE<DFE>(dfe, g);
}

// Turn a DFE into a Sergey DFE, i.e., s(X) = s(0)*(1-g*X)
template <class DFE>
class SergeyDFE{
    public:
        DFE dfe; 
        const double U0;       
        const double g; 
        SergeyDFE(DFE dfe, double g): dfe(dfe), U0(dfe.U0), g(g) {};      
        double get_mutation_rate(Clone const & individual) const { return dfe.get_mutation_rate(individual); }; 
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return std::pow(dfe.get_fitness_effect(random, individual), 1-std::log(individual.fitness)*g );
        };
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_sergey_dfe(DFE dfe, double g) -> SergeyDFE<DFE> {
    return SergeyDFE<DFE>(dfe, g);
}

// Turn a DFE into an Elizabeth DFE, i.e., s(X) = s(0)*exp(-g*X), U(X) = U(0)*exp(g*X)
template <class DFE>
class ElizabethDFE{
    public:
        DFE dfe;
        const double U0;
        const double g; 
        const double gm;
        ElizabethDFE(DFE dfe, double g, double x): dfe(dfe), U0(dfe.U0), g(g), gm(g*std::exp(-x)/(1+std::exp(-x))) {};       
        double get_fitness_effect(Random & random, Clone const & individual) const {
            return std::pow(dfe.get_fitness_effect(random, individual), std::pow(individual.fitness,-g));
        };
        double get_mutation_rate(Clone & individual) const { 
            if(individual.epistasis_parameter < 0){
                individual.epistasis_parameter = dfe.get_mutation_rate(individual)*std::pow(individual.fitness,gm); 
                if(individual.epistasis_parameter > UMAX){
                    individual.epistasis_parameter = UMAX;
                }
            }
            return individual.epistasis_parameter;
        };
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_elizabeth_dfe(DFE dfe, double g, double x) -> ElizabethDFE<DFE> {
    return ElizabethDFE<DFE>(dfe, g, x);
}


// Wrapper class so that the identities of individual mutations can be tracked
template <class DFE>
class TrackedDFE{
    public:
        DFE dfe;   
        const double U0;   
        LabelGenerator label_generator;
        TrackedDFE(DFE dfe): dfe(dfe), U0(dfe.U0) {};       
        double get_mutation_rate(Clone & individual) const { return dfe.get_mutation_rate(individual); };
        Mutation get_mutation(Random & random, Clone & individual, int t=-1) { 
            return Mutation{label_generator.get_next_label(), dfe.get_fitness_effect(random, individual), t};
        };
        void fix_mutation(Mutation const & mutation){};
};

// A utility function for generating a diminishing DFE with minimal syntactic overhead
template <class DFE>
inline auto create_tracked_dfe(DFE dfe) -> TrackedDFE<DFE> {
    return TrackedDFE<DFE>(dfe);
}

// Turn a DFE into a finite sites DFE
template <class DFE>
class FiniteSitesDFE{
    public:
        DFE dfe;
        const int L;
        const double U0;
        std::vector<double> fitness_effects;
        std::vector<bool> sign_flips;
        std::uniform_int_distribution<> draw_label;
        
        FiniteSitesDFE(Random & random, DFE dfe, int L): dfe(dfe), U0(dfe.U0), L(L), draw_label(0,L-1), fitness_effects(L,1), sign_flips(L,false){
        
            for(auto & W : fitness_effects){
                W = dfe.get_fitness_effect(random, Clone());
            }
        };
        
        void fix_mutation(Mutation const & mutation){
            if(mutation.label >= 0 && mutation.label < L){

                sign_flips[mutation.label] = !sign_flips[mutation.label];
            }
        }
        
        double get_mutation_rate(Clone & individual) const { return dfe.get_mutation_rate(individual); };
        
        Mutation get_mutation(Random & random, Clone & individual, int t=-1) { 
                
            auto label = draw_label(random);
            
            bool sign_flip = sign_flips[label];
            for(Mutation & mutation : individual.mutations){
                if(mutation.label == label){
                    sign_flip = !sign_flip;
                }
            }
            
            if(sign_flip)
                return Mutation{label, 1.0/fitness_effects[label], t};
            else
                return Mutation{label, fitness_effects[label], t};
        };
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_finite_sites_dfe(Random & random, DFE dfe, int L) -> FiniteSitesDFE<DFE> {
    return FiniteSitesDFE<DFE>(random, dfe, L);
}

// A beta-parameterized DFE with the density function
//
// rho(x) ~ exp( - (x/s)^beta )
//
// For beta=1.0, we recover the exponential distribution
/*class BetaTrackedDFE{
    public:
       BetaTrackedDFE(double N, double U, double s, double beta): U(U), s(s), beta(beta), mutated(create_bernoulli(U)), draw_num_mutations(create_poisson(N*U)), label_generator_ptr(std::make_shared<LabelGenerator>()) {
           gamma_k = 1.0/beta;
           gamma_theta = std::pow(s,beta);
           draw_s_tothe_beta = create_gamma(gamma_k, gamma_theta);
       };

       decltype(create_bernoulli(0)) mutated;
       decltype(create_poisson(0)) draw_num_mutations;
       decltype(create_gamma(0,0)) draw_s_tothe_beta;

       double get_fitness_effect(Random & random, Clone const & individual) { 
           return std::exp(std::pow(draw_s_tothe_beta(random), gamma_k));
       };

       Mutation get_mutation(Random & random, Clone const & individual) { 
           return Mutation{label_generator_ptr->get_next_label(), get_fitness_effect(random,individual)};
       };
       double U;
       double s;
       double beta;
       double gamma_k;
       double gamma_theta;
       std::shared_ptr<LabelGenerator> label_generator_ptr;  
};*/    


#endif
