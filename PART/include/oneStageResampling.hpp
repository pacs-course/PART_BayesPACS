#ifndef H_ONESTAGERESAMPLING_HPP
#define H_ONESTAGERESAMPLING_HPP

#include <armadillo>
#include <vector>
using namespace arma;

struct Parameters;

/*! @file oneStageResampling.cpp
@brief Performs partitioning and one-stage aggregation step to determine the aggregated posterior's estimate and then draw MCMC samples from it.

@param subchains (\c vector \c < \c const \c mat \c * \c >): vector of pointers to matrices of MCMC samples.
@param par (\c Parameters): Parameters object, storing information on input parameters.
@param min_num_points_block (\c uword): minimum number of MCMC samples falling into each partition's block.
@param y (\c mat): (n_samples x d)-matrix storing MCMC samples drawn from the aggregated posterior.

@par Partitioning phase
After building an object of class RawMCMC based on MCMC samples contained in \c subchains, a \c forest vector representing possible partitions of the sampling area is grown.\n
Each tree in the forest is represented as a vector of leaves. Each leaf node represents an element of the partition of the sampling area.

@par Aggregation and resampling phase
The function calls posteriorResampling() to draw samples from the aggregated posterior, based on estimates contained in the aforementioned forest.

@see Parameters for additional details on input parameters.
@see RawMCMC for additional details on how to handle MCMC samples.
@see Tree for additional details on how to grow a tree in which each leaf node represents an element of the sampling area's partition.
@see MCestimate for additional details on how to handle statistics regarding MCMC samples.
@see posteriorResampling() for additional information on how to draw samples from the aggregated posterior.
*/
void oneStageResampling(const std::vector<const mat*> & subchains, const Parameters & par, const uword min_num_points_block, mat & y);
#endif