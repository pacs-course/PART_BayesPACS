#ifndef H_MATCHSUBSETS_HPP
#define H_MATCHSUBSETS_HPP

#include <armadillo>
#include <vector>
#include <unordered_map>

using namespace arma;

/*! @file matchSubsets.cpp
@brief Determine matching among identifiers of available MCMC samples' matrices (according to \c improve_matching parameter).

@param improve_matching (\c bool): \c true if subsets' indices need to be matched according to median value of squared distance between mean of MCMC samples.
@param MCdraws (\c unordered_map \c < \c uword \c , \c mat \c >): \a key = subset index; \a value = matrix of MCMC samples associated to subset index in \a key.
@param match (\c vector \c < \c vector \c < \c uword \c > \c >): floor(MCdraws_left.size()/2.0) components; match[g] stores indices of subsets belonging to the same group \c g (for all g = 0: floor(MCdraws_left.size()/2.0)-1 ); OBSERVE: if the number of subsets to be matched is odd, the last group will contain three indices of subsets instead of just a couple.

This function is called only if \c pairwise_aggregation is \c true and number of original subsets is at least four.
Details on \c improve_matching parameter:
- If \c improve_matching is \c false, the function matches subsets' indices in pairs without particular ordering. If input number of subsets is odd, last group will contain three elements instead of two.
- If \c improve_matching is \c true, the function matches subsets according to the following metric:
	for each \c s, fixed subset's index that has not yet been matched:
	- compute the following quantity for each \c n subset's index that has not yet been matched: sum_dis_sq(n) = sum(pow(mean(MCdraws[n],0) - mean(MCdraws[s],0), 2),1)
	- match subset with index \c s to subset with index \c n such that sum_dis_sq(n) = median(sum_dis_sq)
*/
void matchSubsets(const bool improve_matching, const std::unordered_map<uword,mat> & MCdraws, const bool verbose, std::vector<std::vector<uword> > & match);
#endif