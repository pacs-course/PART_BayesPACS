#include <iostream>
#include <cassert>
#include "Parameters.hpp"
#include "RawMCMC.hpp"
#include "MCestimate.hpp"
#include "Tree.hpp"
#include "posteriorResampling.hpp"

using namespace arma;

void oneStageResampling(const std::vector<const mat*> & subchains, const Parameters & par, const uword min_num_points_block, mat & y){

	std::cout << "******************* Grow a forest of " << par.ntree << " trees *******************" << std::endl << std::endl;
	
	const RawMCMC rawMCMC(subchains); // store information on MCMC samples
	
	const uword nMCsamples = rawMCMC.get_samples().n_rows;
	assert(nMCsamples>0 && "Number of samples must be at least 1!!");
	
	// OBSERVE: If nMCsamples<=min_num_points_block each tree in the forest will have only the root node, storing information on rawMCMC object.
	
	const uvec dim_cut = find(rawMCMC.get_area().row(1)-rawMCMC.get_area().row(0)>par.min_cut_length); // dimensions along which cutting point should be searched are those that do not violate minimum length of area
	
	if(par.verbose){
		std::cout << "Number of samples: " << nMCsamples << std::endl;
		rawMCMC.get_area().print("Sampling Area: ");
		dim_cut.print("Dimensions along which cutting point will be searched: ");
		std::cout << "\n Minimum number of samples in each partition's block: " << min_num_points_block << std::endl;
	}
	
	std::vector<std::vector<MCestimate> > forest(par.ntree);
	for(uword n=0; n<par.ntree; n++){ // forest[n] = vector of trees; each tree is represented as vector of leaf nodes; each leaf represents an element of the partition of the sampling area		
		
		// Create Tree object for space partitioning
		Tree tree(par.kd_cut, min_num_points_block, par.min_cut_length);
		tree.grow(rawMCMC,dim_cut, par.verbose, forest[n]); // grow random partition tree
		
		if(par.verbose){
			std::cout << "Dimension of " << n <<"-th tree : " << forest[n].size() << std::endl << std::endl;	
		}
		
	}
	
	std::cout << "******************* Draw MCMC samples from the aggregated posterior *******************" << std::endl << std::endl;			
	posteriorResampling(forest, par.n_samples, par.gaussian_smooth, par.verbose, y);

};