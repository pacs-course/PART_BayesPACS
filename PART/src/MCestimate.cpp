#include <stdexcept>
#include <cassert>
#include <cmath>
#include "MCestimate.hpp"
#include "RawMCMC.hpp"
#include "computeMeanCov.hpp"

using namespace arma;

// Constructor
MCestimate::MCestimate(const RawMCMC & rawMCMC, const uvec & idx, const mat & area):sampling_area(area){

	const uword nMCsamples = idx.n_rows;
	assert(nMCsamples>0 && "Number of samples must be at least 1!!");

	const uword M = rawMCMC.get_n_subsets(); // Number of subsets
	assert(M>0 && "Number of subsets cannot be less than 1!!!");
	
	const uvec & mark = rawMCMC.get_mark()(idx);
	
	const double sum_logN = rawMCMC.get_sumLogN(); // = \log( \prod_{m=0}^{M-1} N(m) ) = sum_{m=0}^{M-1} \log(N(m))
	if(M>1){ // unnormalized density at current node: prob = \prod_{m=0}^{M-1} \{ \frac{nBlock(m)}{N(m)} \} / (|area|^{M-1})
		
		double sum_log_nBlock=0.0; // = \log( \prod_{m=0}^{M-1} nBlock(m) ) = \sum_{m=0}^{M-1} \log( nBlock(m) ) )
		for(uword m=0; m<M; m++){
			const uvec id_sub = find(mark==m); // id_sub.n_rows: number of samples belonging to subset m falling into current node
			sum_log_nBlock += std::log(static_cast<double>(id_sub.n_rows+0.1)); 
		}
		
		double sum_log=0.0; // \log(|area|^{M-1}) = (M-1) * \log(|area|) = (M-1) * \log( \prod_{p=0}^{d-1}( area(1,p)-area(0,p) ) ) = (M-1) * \sum_{p=0}^{d-1}( \log(area(1,p)-area(0,p)) )
		const rowvec & length_area = area.row(1) - area.row(0); // length of area along each direction p=0:d-1
		for(const double l: length_area){
			sum_log += std::log(l);
		}
		
		log_prob = sum_log_nBlock - sum_logN - (M-1)*sum_log;
		
	}else{ // unnormalized density at current node: prob = nBlock/N = nMCsamples/N
		log_prob = std::log(nMCsamples) - sum_logN;
	}
	
	// compute mean and covariance of MCMC samples in current node
	computeMeanCov(rawMCMC.get_samples().rows(idx), mark, M, mean, Cov);
	
};