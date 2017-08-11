#ifndef H_PAIRWISEAGGREGATION_HPP
#define H_PAIRWISEAGGREGATION_HPP

#include <armadillo>
#include <vector>
using namespace arma;

struct Parameters;

/*! @file pairwiseAggregation.cpp
@brief Performs the aggregation step combining MCMC samples' matrices in groups of two or three and determining, for each of these groups, an estimate of the aggregated posterior.

@param par (\c Parameters): \c Parameters object, storing information on input parameters.
@param subchains (\c vector \c < \c mat \c >): vector of \c par.M matrices of MCMC samples; each m_th element is a matrix of dimensions (\c par.N[m] x \c d); (\c d = number of attributes).
@param aggr_post_samples (\c mat): (\c par.n_samples x \c d)-matrix storing MCMC samples drawn from the aggregated posterior.

This function is called when \c pairwise_aggregation is \c true and \c par.M > 3.

The function repeates the following procedure for  floor(log(par.M)/log(2)) stages:
- calls matchSubsets() to determine an appropriate matching among MCMC samples' matrices. The matching's result is a list of groups, each of which contains two (or at most three) indices of matrices.
- calls oneStageResampling() on each group of two (or at most three) matrices. The resulting estimated posterior for each group is used as input for next stage.

The whole procedure terminates with a common aggregated posterior.

@see Parameters for additional details on input parameters.
@see matchSubsets() for additional details on how to determine an appropriate matching among subsets.
@see oneStageResampling() for additional details on how to determine an estimate of the aggregated posterior and resample from it.
*/
void pairwiseAggregation(const Parameters & par, const std::vector<mat> & subchains, mat & aggr_post_samples);
#endif