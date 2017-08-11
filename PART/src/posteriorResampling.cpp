#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include "MCestimate.hpp"
#include "posteriorResampling.hpp"

using namespace arma;

void posteriorResampling(const std::vector<std::vector<MCestimate> > & forest, const uword n_samples, const bool gaussian_smooth, const bool verbose, mat & y){
	
	const uword ntree = static_cast<uword>(forest.size()); // number of trees in forest
	
	// uniformly sample tree indices for n_samples-times
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> unif_int(0,ntree-1);
	uvec sampled_tree(n_samples);
	for(uword n=0; n<n_samples; n++){
		sampled_tree(n) = static_cast<uword>(unif_int(gen));
	}
	
	// count how many times each tree has been selected
	uvec tree_counted(ntree);
	for(uword n=0; n<ntree; n++){
		const uvec tmp = find(sampled_tree==n);
		tree_counted(n) = tmp.n_rows;
	}
	
	const uword d = forest[0][0].mean.n_cols;
	if(y.n_rows != n_samples || y.n_cols != d){
		y.set_size(n_samples,d);
	}
	
	uword num_sampled = 0; // number of MCMC samples drawn from the aggregated posterior
	for(unsigned int ctree=0; ctree<static_cast<unsigned int>(ntree); ctree++){
	
		if(verbose){
			std::cout << "Index of selected tree :" << ctree << std::endl << std::endl;
		}

		// Select the current tree from which resampling of leaf nodes should be done
        const std::vector<MCestimate> tree(forest[ctree]);
		
        /* Sample (with replacement) leaf nodes in current tree for (tree_counted(ctree))-times with weigths == prob */
		std::vector<double> prob(tree.size()); // weight of each mixture component
		unsigned int pos=0;
		for(auto const & i: tree){
			prob[pos] = std::exp(i.log_prob);
			pos++;
		}
		
		const uword freq_tree = tree_counted(ctree); // how many times current tree have been sampled
		
		uvec sampled_nodes(freq_tree);
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<> discr(prob.begin(), prob.end());
		for(uword n=0; n<freq_tree; n++){
			sampled_nodes(n) = static_cast<uword>(discr(gen));
		}
		
		const uvec unique_nodes = unique(sampled_nodes); // select unique values of sampled_nodes in ascending order
		const uword n_unique = unique_nodes.n_rows; // number of unique sampled_nodes
		
		uvec nodes_count(n_unique);// how many times each node in current tree had been sampled
		pos=0;
		for(const uword i: unique_nodes){
			const uvec tmp = find(sampled_nodes==i);
			nodes_count(pos) = tmp.n_rows;
			pos++;
		}
		
		if(verbose){
			std::cout << "Indices of nodes onto which resampling is performed: ";
		}
		
        for(uword node_index = 0; node_index<n_unique; node_index++){ // for each unique value of sampled nodes
			
			// select the leaf node onto which resampling should be performed
			const unsigned int id = static_cast<unsigned int>(unique_nodes(node_index));
			const MCestimate leaf = tree[id]; 

			if(verbose){
				std::cout << id << "    ";
			}
			
			const uword freq_node = nodes_count(node_index); // how many times current node have been sampled in current tree
			
			vec eigval;
			mat eigvec;
			if(gaussian_smooth && (leaf.Cov).is_finite() && eig_sym(eigval, eigvec, leaf.Cov) && eigval.is_finite() && eigvec.is_finite() && all(eigval>0.0) ){ // covariance positive definite and invertible
				// Generate freq_node samples from Gauss(leaf.mean,leaf.cov)
				
				eigval = pow(eigval,0.5);
				assert(eigval.is_finite() && "Error in sqrt of eigenval!!");
				
				y.rows(num_sampled,num_sampled+freq_node-1) = randn<mat>(freq_node,d)*diagmat(eigval)*(eigvec.t()) + repmat(leaf.mean,freq_node,1);
				
			}else{ 
				// Generate freq_node samples from Unif(leaf.sampling_area)
				
				y.rows(num_sampled,num_sampled+freq_node-1)=repmat(leaf.sampling_area.row(1)-leaf.sampling_area.row(0),freq_node,1)%randu<mat>(freq_node,d) + repmat(leaf.sampling_area.row(0),freq_node,1);
				
			}
			
			num_sampled += freq_node;// update number of samples drawn from the estimated posterior

		}
		
		if(verbose){
			std::cout << std::endl;
		}
		
	}
	std::cout << "Posterior resampling DONE! " << std::endl << std::endl;
	
};