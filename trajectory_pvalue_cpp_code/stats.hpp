#ifndef STATS_HPP
#define STATS_HPP

#include <random>
#include <functional>
#include <iostream>

typedef std::mt19937 Random;
static std::uniform_real_distribution<> sample_uniform(0,1);
static std::normal_distribution<> sample_normal(0,1);
static std::exponential_distribution<> sample_exponential(1);

inline Random create_random(unsigned int seed=0, bool print_seed=false){
    if(seed==0){
        std::random_device rd;
        seed = rd();
    }
    if(print_seed){
        std::cout << seed << std::endl;
    }
    return Random(seed);
}

inline std::uniform_int_distribution<> create_random_int(int min, int max){ return std::uniform_int_distribution<>(min,max);} 
inline std::uniform_int_distribution<> create_uniform_int(int min, int max){ return std::uniform_int_distribution<>(min,max);} 
inline std::poisson_distribution<> create_poisson(double lam) { return std::poisson_distribution<>(lam);}
inline std::bernoulli_distribution create_bernoulli(double p) { return std::bernoulli_distribution(p); }
inline std::gamma_distribution<> create_gamma(double k, double theta){ return std::gamma_distribution<>(k,theta); }; /* x^(k-1)*exp(-x/theta) */

inline int sample_precomputed_small_poisson(Random & random, double L){
    double p = sample_uniform(random);
    int k=0;
    while(p > L){
        ++k;
        p*=sample_uniform(random);
    }
    return k;
}

inline int sample_small_poisson(Random & random, double mean){
    double L = exp(-mean);
    return sample_precomputed_small_poisson(random,L);
}

/*function gammaln for log of gamma function*/
inline double gammaln(double xx) {
    double x,y,tmp,ser;
    double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    long j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
} /*gammaln*/


inline int sample_poisson(Random & random, double lam){
    if(lam >= 10.0){
        int k;
        double U, V, slam, loglam, a, b, invalpha, vr, us;  
        slam = sqrt(lam);
        loglam = log(lam);
        b = 0.931 + 2.53*slam;
        a = -0.059 + 0.02483*b;
        invalpha = 1.1239 + 1.1328/(b-3.4);
        vr = 0.9277 - 3.6224/(b-2); 
        while (1){
            U = sample_uniform(random) - 0.5;
            V = sample_uniform(random);
            us = 0.5 - fabs(U);
                k = (int) floor((2*a/us + b)*U + lam + 0.43);
                if ((us >= 0.07) && (V <= vr)){
                    return k;
                }
                if ((k < 0) || ((us < 0.013) && (V > us))){
                    continue;
                }
                if ((log(V) + log(invalpha) - log(a/(us*us)+b)) <= (-lam + k*loglam - gammaln(k+1))){
                    return k;
                }      
        }   
    }
    else{
        return sample_small_poisson(random, lam);
    }
}

#endif
