#ifndef H_PARAMETERS_HPP
#define H_PARAMETERS_HPP

#include <armadillo>
#include <string>
using namespace arma;

/*! @file Parameters.cpp
@brief Basic class to handle information regarding input parameters.
*/
struct Parameters{
	
	uword M; /*!< Number of subsets.*/
	bool pairwise_aggregation; /*!< \c true if pairwise aggregation method to be used. */
	bool improve_matching; /*!< \c true if in pairwise aggregation subsets' indices need to be matched according to median value of squared distance between mean of MCMC samples. OBSERVE: \c improve_matching is ignored if \c pairwise_aggregation is \c false. */
	bool halve; /*!< \c true if \c min_frac_points_block need to be halved at each step of pairwise aggregation. OBSERVE: \c halve is ignored if \c pairwise_aggregation is \c false. */
	uword ntree; /*!< Number of trees to be grown in forest. */
	bool kd_cut; /*!< \c true if cutting point determined as median of MCMC samples; \c false if cutting point maximizes the empirical log-likelihood. */
	double min_frac_points_block; /*!< Minimum fraction of MCMC samples falling into each partition's block. */
	double min_cut_length; /*!< Minimum length of each dimension of partition's block. */
	uword resample_N; /*!< Minimum number of MCMC samples to be drawn from the aggregated posterior at each step. */
	uword n_samples; /*!< Number of MCMC samples to be drawn from the aggregated posterior. */	
	bool gaussian_smooth; /*!< \c true if the aggregated posterior density is a mixture of local gaussian densities. */
	bool verbose; /*!< \c true if all messages must be displayed; =0 if only basic information displayed. */	
	uword d; /*!< Number of attributes. */
	uvec N;	/*!< M-vector such that N(m) = number of MCMC samples stored on subset m. */

	//! Constructor from parameters' file and number of attributes in current problem
	Parameters(const std::string & paramsFile);
	
	//! Display information on currently used parameters on screen
	void display() const;
	
};
#endif