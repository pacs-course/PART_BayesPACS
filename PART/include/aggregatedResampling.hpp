#ifndef H_AGGREGATEDRESAMPLING_HPP
#define H_AGGREGATEDRESAMPLING_HPP

#include <armadillo>
#include <vector>
using namespace arma;

struct Parameters;

/*! @file aggregatedResampling.cpp
@brief Defines how to perform the aggregation step (pairwise or one-stage) to determine an estimate of the posterior distribution and to draw samples from it.

@param par (\c Parameters): \c Parameters object, storing information on PART parameters.
@param subchains (\c vector \c < \c mat \c >): vector of \c par.M matrices of MCMC samples; each m_th element is a matrix of dimensions (\c par.N(m) x \c d); (\c d := number of features).
@param aggr_post_samples (\c mat): (\c par.n_samples x \c d)-matrix storing MCMC samples drawn from the aggregated posterior.

@par Details on \c pairwise_aggregation parameter:
If \c pairwise_aggregation is \c false, the posterior estimate is obtained combining all MCMC samples' matrices together in one single step and drawing samples from the resulting aggregated posterior. This is performed by oneStageResampling().
\n
If \c pairwise_aggregation is \c true, the posterior estimate is obtained calling function pairwiseAggregation() that combines the input MCMC samples in groups of two or three subsets.
\n The matching among subsets is determined through matchSubsets().
\n OBSERVE: if \c par.M is not larger than three, even if \c pairwise_aggregation is \c true, all subsets need to be matched together hence the function oneStageResampling() is called directly.

@see Parameters for additional details on input parameters.
@see oneStageResampling() for additional details on how to determine an estimate of the aggregated posterior and resample from it in one single stage.
@see pairwiseAggregation() for additional details on how to determine an estimate of the aggregated posterior and resample from it using multiple stages.
@see matchSubsets() for additional details on how to determine an appropriate matching among subsets.
*/
void aggregatedResampling(const Parameters & par, const std::vector<mat> & subchains, mat & aggr_post_samples);
#endif