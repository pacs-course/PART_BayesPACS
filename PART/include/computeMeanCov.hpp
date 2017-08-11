#ifndef H_COMPUTEMEANCOV_HPP
#define H_COMPUTEMEANCOV_HPP

#include <armadillo>
using namespace arma;

//! Compute mean and covariance of MCMC samples
/*!
If samples belong to different subsets, this function computes the mean and covariance using block estimates.\n
If the estimate of the covariance matrix obtained through block estimates is not positive definite or the number of subsets in current object is exactly one, simple estimates are used.
*/
void computeMeanCov (const mat & samples, const uvec & mark, const uword M, rowvec & mu, mat & Sigma);
#endif