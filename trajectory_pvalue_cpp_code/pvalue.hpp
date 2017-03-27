#ifndef PVALUE_HPP
#define PVALUE_HPP

#include <cmath>

// Helper functions for calculating empirical pvalues

const int min_numerator_counts = 100;
const double pvalue_threshold = std::pow(0.05,1.0/3);

inline double clip_pvalue(double pvalue){
    if(pvalue>pvalue_threshold){
        return 1.0;
    }
    else{
        return pvalue;
    }
}

inline double calculate_combined_pvalue(const double p1, const double p2, const double p3){
    return clip_pvalue(p1)*clip_pvalue(p2)*clip_pvalue(p3);
}

inline double calculate_pvalue_from_counts(int num_greater, int total_num){
    return (num_greater+1.0)/(total_num+1.0);
}

template<class T1, class T2> inline double calculate_pvalue_from_sorted_list(T1 const & sorted_list, T2 const & query_val){

    int num_greater = sorted_list.end()-std::lower_bound(sorted_list.begin(), sorted_list.end(), query_val);
    int total_num = sorted_list.size();
    
    return calculate_pvalue_from_counts(num_greater, total_num);
}



#endif