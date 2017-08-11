#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "Parameters.hpp"
#include "aggregatedResampling.hpp"
#include "pairwiseAggregation.hpp"
#include "oneStageResampling.hpp"

using namespace arma;

void aggregatedResampling(const Parameters & par, const std::vector<mat> & subchains, mat & aggr_post_samples){
	
	assert(par.M==subchains.size() && "Number of subsets stored in input matrix and input number of subsets must coincide!!");
	
	if(par.pairwise_aggregation && par.M>3){
		
		std::cout << "******************* Start pairwise aggregation phase *******************" << std::endl << std::endl;
		pairwiseAggregation(par, subchains, aggr_post_samples);
		
	}else{
		// Even if par.pairwise_aggregation==true, if par.M<=3 no need to search for subsets' matching.
		// Aggregate samples from all subsets to perform partitioning and resampling.
		
		std::cout << "******************* Start one-stage aggregation phase *******************" << std::endl << std::endl;

		std::vector<const mat*> ptrSubchains(par.M,nullptr); // vector of pointers to subchains content
		for(unsigned int m=0; m<par.M; m++){
			ptrSubchains[m] = &(subchains[m]);
		}
		oneStageResampling(ptrSubchains, par, static_cast<uword>(std::ceil(max(par.N)*par.min_frac_points_block)), aggr_post_samples);
	}
	
};

