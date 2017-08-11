#ifndef H_MCESTIMATE_HPP
#define H_MCESTIMATE_HPP

#include <armadillo>

using namespace arma;

class RawMCMC;

/*! @file MCestimate.cpp
@brief Basic struct to handle information regarding each leaf node, representing one element of the partition of sampling area.

For each partition's block, we store information on:
- \c sampling_area at current leaf,
- \c mean, \c Cov: mean and covariance of MCMC samples falling into current partition's block,
- \c log_prob: logarithmic form of the unnormalized estimate of aggregated posterior's component associated to current partition's block.\n
	In logarithmic form, the unnormalized estimate of the aggregated posterior is:\n
	\f[ \log(prob) = sum(\log(n_{samples})) - sum(\log(N)) - (M-1)*sum(\log(area.row(1)-area.row(0))) \f]
	\n where:
	- N (\c uvec) such that N(m) = number of MCMC samples stored on machine m.
	- sampling_area (\c mat) coincide with sampling area in current partition's block.
	- n_samples (\c uvec) such that n_samples(m) = number of MCMC samples in current block that have been drawn on machine m.

@see RawMCMC for additional details on how to store MCMC samples.
@see Tree, oneStageResampling() and posteriorResampling() for additional details on how to use information stored in MCestimate.
*/
struct MCestimate{
	
	double log_prob; //!< Unnormalized estimate of aggregated posterior's component associated to current partition's block.
	mat sampling_area; //!< Sampling area at current partition's block.
	rowvec mean; //!< Mean of samples falling into current partition's block.
	mat Cov; //!< Covariance of samples falling into current partition's block.
	
	//! Constructor
	MCestimate(const RawMCMC & rawMCMC, const uvec & idx, const mat & area);
	
};
#endif