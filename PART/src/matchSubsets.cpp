#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <list>
#include <iterator>
#include "matchSubsets.hpp"

using namespace arma;

void matchSubsets(const bool improve_matching, const std::unordered_map<uword,mat> & MCdraws, const bool verbose, std::vector<std::vector<uword> > & match){
	
	const uword subs_left = MCdraws.size(); // Number of subsets to be matched
	assert(subs_left>3 && "matchSubsets called also when MCdraws.size()<=3!!!");
	
	const uword G = static_cast<uword>(std::floor(subs_left/2.0)); // Number of groups after matching		
	match.resize(G, std::vector<uword>(2)); // match[g]: subsets' indices that belong to the same group g
	
	// Randomness in matching subsets: randomly shuffle keys stored in MCdraws
	uvec id_sub(subs_left);
	uword pos=0;
	for(auto const & it: MCdraws){
		id_sub(pos)=it.first;
		pos++;
	}
	id_sub = shuffle(id_sub);// randomly shuffle
	
	if(improve_matching){ // match subsets according to the following metric:
		/* for each fixed subset's index s that has not yet been matched:
			- compute:	sum_dis_sq[n] = accu(pow(mean(MCdraws[n],0) - mean(MCdraws[s],0),2))	for all subsets' index n that has not yet been matched
				where
					mean(matrix,0) = mean of matrix on row indices
					pow(vec,2) = element-wise power to 2 of vec components
					accu(vec) = sum all components of vec
			- match subset's index s to subset's index n corresponding to median(sum_dis_sq), that is:	sum_dis_sq[n] = sum_dis_sq(sort_index(sum_dis_sq)(static_cast<uword>(std::floor(sum_dis_sq.n_rows/2.0))));
		*/
		
		if(verbose){
			std::cout << "Use improved matching: associate subsets according to median sum of squared distance among mean of MCMC samples ... " << std::endl << std::endl;
		}
		
		std::list<uword> subs_left_set; // stores indices of subsets that still need to be matched
		for(const uword s: id_sub){
			subs_left_set.emplace_back(s);
		}
		
		uword n_groups = 0; // number of groups formed
		while(subs_left_set.size()>3){ // until there are no more than 3 subsets left out from matching
			
			const uword keyF = *subs_left_set.begin(); // key associated to element that need to be matched
			
			auto const itF = MCdraws.find(keyF); // iterator to MCdraws' element corresponding to keyF
			assert(itF != MCdraws.end() && "ERROR: keyF does not exist in MCdraws!!");

			const rowvec meanF = mean(itF->second,0); // mean (by row) of matrix associated to keyF

			vec sum_dis_sq(subs_left_set.size()); // store sum of distances squared for each element in subs_left_set
			uword pos=0;
			for(const uword s: subs_left_set){
				
				auto const it = MCdraws.find(s);
				assert(it != MCdraws.end() && "ERROR: key does not exist in MCdraws!!");

				const rowvec dis(mean(it->second,0)-meanF); // mean distance w.r.t. matrix associated to keyF
				
				sum_dis_sq(pos) = accu(dis%dis); // sum of pow 2 of each element of dis
				pos++;
			}

			const uvec id_sorted = sort_index(sum_dis_sq); // id of sum_dis_sq sorted from lowest to highest value
			
			auto itS = subs_left_set.begin();
			std::advance(itS, static_cast<int>(id_sorted(static_cast<uword>(std::floor(subs_left_set.size()/2.0))))); // iterator to second index in current group
			
			// insert pair into match vector
			match[n_groups][0] = keyF;
			match[n_groups][1] = *itS;
			n_groups++;
			
			// remove form subs_left_set id of subsets that have already been matched
			subs_left_set.erase(itS);
			subs_left_set.erase(subs_left_set.begin());
			
		}
		// At this point we have subs_left_set.size() <= 3, meaning that we need to create the last group with all subsets currently stored in subs_left_set
		assert(G-1==n_groups && "!!! ERROR: number of pairs is different from number of groups in matching !!!");
		match[n_groups].resize(subs_left_set.size());
		std::copy(subs_left_set.begin(),subs_left_set.end(), match[n_groups].begin());
		
	}else{ // match subsets without specific order

		if(verbose){
			std::cout << "Match subsets without specific order ... " << std::endl << std::endl;	
		}
		
		// match pairs of indices
		for(uword g=0; g<G; g++){
			match[g][0]=id_sub(2*g);
			match[g][1]=id_sub(2*g+1);
		}
		
		// if the number of input subsets is odd, add the index of subset still left out to the last group
		if(subs_left>2*G){
			match.back().push_back(id_sub(subs_left-1));
		}
		
	}
	
	if(verbose){
		std::cout << " ---> Resulting matching: " << std::endl;
		for(const auto & g: match){
			std::cout << "		";
			for(const uword s: g){
				std::cout << s << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
};