#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include "pairwiseAggregation.hpp"
#include "Parameters.hpp"
#include "matchSubsets.hpp"
#include "oneStageResampling.hpp"

using namespace arma;

void pairwiseAggregation(const Parameters & par, const std::vector<mat> & subchains, mat & aggr_post_samples){
	
	const uword M = subchains.size(); // Number of subsets
	assert(M>3 && "!!! ERROR: Number of subsets in pairwise aggregation need to be at least 4!!!");
	
	const uword stages = static_cast<uword>(std::floor(std::log(M)/std::log(2)));
	if(par.verbose){
		std::cout << "Number of stages to conclude pairwise aggregation: " << stages << std::endl << std::endl;
	}
	
	double min_frac_points_block = par.min_frac_points_block; // minimum fraction of MCMC samples in each partition's block
	
	// if par.halve==true, need to halve at each stage the fraction of MCMC samples falling into each partition's block 
	if(par.halve && stages>1){// update minimum fraction of points in each block to have at final stage the value par.min_frac_points_block
		min_frac_points_block = static_cast<double>(par.min_frac_points_block * std::pow(2.0,stages-1));
	} 
	
	uword min_num_points_block = static_cast<uword>(std::max(static_cast<int>(std::ceil(max(par.N)*min_frac_points_block)),3));
	assert(min_num_points_block>2 && "At least three samples in each partition's block!!");
	
	if(par.verbose){
		std::cout << "Minimum number of MCMC samples in each partition's block at first aggregation's step: " << min_num_points_block << std::endl << std::endl;
	}
	
	std::cout << "******************* Start matching subsets *******************" << std::endl << std::endl;
	
	std::unordered_map<uword,mat> MCdraws; // initialize MCdraws with input matrices
	for(uword i=0; i<M; i++){
		auto const it = MCdraws.emplace(i,subchains[i]);
		assert(it.second && "Emplace in unordered_map failed!!");
	}
	
	while(MCdraws.size()>1){ // until there is only one matrix left
		
		std::vector<std::vector<uword> > match; // match[g]: indices of subsets that belong to the same group g; match[g] always stores 2 elements, except for last group that could contain at most 3 elements
		if(MCdraws.size()>3){ // if more than three subsets left out, determine matching among subsets' indices
			
			matchSubsets(par.improve_matching, MCdraws, par.verbose, match);
			
		}else{ // the number of subsets left out is not larger than three, hence all subsets will be grouped together
			
			match.resize(1);
			for(auto const & i: MCdraws){
				match.front().push_back(i.first);
			}
			
		}
		
		uword id_group=0; // index of group analyzed, just for screen informations
		for(auto const & g: match){ // for each group in match
			
			std::cout << "******************* One-stage aggregation and resampling on group " << id_group << " *******************"<< std::endl << std::endl;
			
			const unsigned int n_matched = g.size(); // number of subsets belonging to current group; either 2 or 3
			assert(n_matched>1 && n_matched<4 && "Error in matching: wrong number of subsets in match!!");
			
			std::vector<const mat*> inner_subchains(n_matched,nullptr); // inner_subchains[i]: pointer to const matrix representing MCMC samples corresponding to subset i in group g
			unsigned int pos=0;
			for(const uword i: g){ // for each subset index i belonging to current group
				auto const it = MCdraws.find(i);
				assert(it != MCdraws.end() && "ERROR: key does not exist in MCdraws!!");
				inner_subchains[pos] = &(it->second); // store matrix associated to subset i
				pos++;
			}
			
			mat post_samples; // store result of posterior resampling associated to current group
			oneStageResampling(inner_subchains, par, min_num_points_block, post_samples);
			MCdraws.find(g[0])->second = std::move(post_samples); // store posterior resampling result in MCdraws' key corresponding to first element of current group
			
			// remove matrices that are no more useful (those corresponding to keys associated to other elements of g except first)
			for(uword i=1; i<n_matched; i++){
				auto const it = MCdraws.erase(g[i]);
				assert(0 != it && "Element not removed from MCdraws!!");
			}
			
			id_group++;
		}
		assert(match.size()==MCdraws.size() && "Error in updating dimension of MCdraws!!");
		// At this point: MCdraws.size() == match.size()
		
		// if par.halve==1, halve minimum fraction of MCMC samples falling into each partition's block
		if(par.halve){
			min_frac_points_block = min_frac_points_block/2.0;
		}
		
		// Update minimum number of MCMC samples falling into each partition's block for next step
		min_num_points_block = static_cast<uword>(std::max(static_cast<int>(std::ceil(par.n_samples*min_frac_points_block)), 3));
		
	}
	assert(1==MCdraws.size() && "Error in MCdraws dimensions!!");
	aggr_post_samples = MCdraws.begin()->second;

};
