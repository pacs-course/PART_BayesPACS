#ifndef H_ONESTAGERESAMPLING_HPP
#define H_ONESTAGERESAMPLING_HPP

#include <armadillo>
#include <vector>
using namespace arma;

class MCestimate;

/*! @file posteriorResampling.cpp
@brief Draws samples from the aggregated posterior.

@param forest (\c vector \c < \c vector \c < \c MCestimate \c > \c >): information regarding each binary tree, representing sampling area partition.
@param n_samples (\c uword): number of MCMC samples to be drawn from the aggregated posterior.
@param gaussian_smooth (\c bool): \c true if the aggregated posterior density is a mixture of local gaussian densities; \c false if kernels are uniforms.
@param verbose (\c bool): \c true if all messages must be displayed; \c false if only basic information displayed.
@param y (\c mat): (\c n_samples x \c d )-matrix storing MCMC samples drawn from the aggregated posterior.

@par Resampling procedure
For each binary tree stored in \c forest, the function samples leaf nodes using as sampling weights the corresponding unnormalized densities (exp(log_prob)).
The number of times each node need to be sampled is given by the probability of selecting the current tree among all possible trees in the forest.
For each unique value of previously sampled nodes, draw MCMC samples from the aggregated posterior.
The number of MCMC samples to be drawn from current node is given by the probability of selecting current node among all possible nodes in current tree.

The criterion according to which the resampling at each node is performed is established through \c gaussian_smooth:
- If \c gaussian_smooth is \c true, local gaussian densities are used to approximate the posterior distribution.
- If \c gaussian_smooth is \c false, local uniform densities are used to approximate the posterior distribution.

OBSERVE: it is possible to use local gaussians ONLY IF the covariance of MCMC samples stored in selected node is a positive definite matrix.
If this condition is not satisfied, local uniforms will be used instead of gaussians.

@see MCestimate for additional details on how to store statistics regarding MCMC samples.				
@see Tree for additional details on how to grow a partition tree.
@see oneStageResampling() for additional details on how to grow a forest of partition trees.
*/
void posteriorResampling(const std::vector<std::vector<MCestimate> > & forest, const uword n_samples, const bool gaussian_smooth, const bool verbose, mat & y);
#endif