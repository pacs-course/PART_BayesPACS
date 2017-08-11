#include <iostream>
#include <cassert>
#include <random>
#include <cmath>
#include "RawMCMC.hpp"
#include "MCestimate.hpp"
#include "Tree.hpp"
#include "computeMeanCov.hpp"

using namespace arma;

// Constructor of binary tree (represented through \a leaves vector)
Tree::Tree(const bool kdCut, const uword minNumPointsBlock, const double minCutLength):kd_cut(kdCut),min_num_points_block(minNumPointsBlock), min_cut_length(minCutLength){};

// getter of kd_cut attribute
const bool Tree::get_kd_cut() const { return kd_cut; };

// getter of min_num_points_block attribute
const uword Tree::get_min_num_points_block() const { return min_num_points_block; };

// getter of min_cut_length attribute
const double Tree::get_min_cut_length() const { return min_cut_length; };

// setter of kd_cut attribute
void Tree::set_kd_cut(const bool kdCut){ this->kd_cut = kdCut; };
	
// setter of min_num_points_block attribute
void Tree::set_min_num_points_block(const uword minNumPointsBlock){ this->min_num_points_block = minNumPointsBlock; };
	
// setter of min_cut_length attribute
void Tree::set_min_cut_length(const double minCutLength){ this->min_cut_length = minCutLength; };

// Constructor of root node
Tree::Node::Node():parent(nullptr),childLeft(nullptr),childRight(nullptr){};

// Method to define attibutes of child node and to link it to current node
void Tree::Node::generateChild(const uword dim, const double cut, const double min_cut_length, const bool is_left, const uvec & idx, Node & node){
	
	node.parent = this; // input node is child of current node hence it need to store current node as parent
	node.indices = this->indices(idx); // define indices as subset of parent's indices
	mat area_child = this->area;
	if(is_left){ // input node is left child of current node
		area_child(1,dim) = cut; // update right boundary of area with cutting point
		this->childLeft = &node; // current node stores input node as left child
	}else{ // input node is right child of current node
		area_child(0,dim) = cut; // update left boundary of area_right with cutting point
		this->childRight = &node; // current node stores input node as right child
	}			
	node.area = area_child;
	node.dim_cut = find(node.area.row(1) - node.area.row(0) > min_cut_length);// pre-processing dimensions to cut for child
	
};

// Method to search for a partition of current area.
void Tree::findCut(const RawMCMC & rawMCMC, const Node & node, double & cut, uword & dim, vec & samples_dim, bool & valid_cut){
	
	valid_cut = false; // true if and only if valid cutting point found
	
	const uword nMCsamples = node.indices.n_rows;
	const mat & samples = rawMCMC.get_samples().rows(node.indices);
	
	// Random devices to select randomly dimension along which cut should be searched
	std::random_device rd;
	std::mt19937 gen(rd()); 
	
	// use this to select candidate dimensions and to remove already checked dimensions
	std::vector<uword> candidate_dims( conv_to< std::vector<uword> >::from(node.dim_cut) );
	while( (candidate_dims.size()>0) && (!valid_cut) ){ // proceed until all dimensions have been checked or a valid partition has been found
		
		// randomly select one of the components of dim_cut that has not yet been analyzed
		std::uniform_int_distribution<> dis(0, candidate_dims.size()-1);
		const uword id_dim = dis(gen);
		dim = candidate_dims[id_dim];

		samples_dim = samples.col(dim); // select only the components along dimension dim of the current MCMC samples
		assert(min(samples_dim)>=node.area(0,dim) && max(samples_dim)<=node.area(1,dim) && "Error in searching cutting point: MCMC samples do not fall in sampling area!!");
		
		// select candidate cutting point according to cutting rule
		if(nMCsamples==1){
			cut = samples_dim(0); // prevent length error in ml_cutting_rule
		}else{
	
			if(kd_cut){ // compute cutting point using kd rule (i.e. cutting point is median point of samples)
				cut = median(samples_dim);
			}else{ // compute cutting point using ML rule
				cut = ml_cutting_rule(rawMCMC,node,dim);
			}
		}
		
		// check validity of cutting point
		valid_cut = checkCut(cut,samples_dim,node.area(0,dim),node.area(1,dim));
		
		candidate_dims.erase(candidate_dims.begin()+static_cast<int>(id_dim)); // remove dimension already checked
		
	}
	// At this point, either all dimensions in dim_cut have been checked or a valid partition has been found
		
	// Add randomness to build different trees
	/* If the following conditions are verified:
	 - node is root
	 - there is only one dimension along which cutting point could be searched
	 - a valid cutting point has already been selected through previous search
	 each partition tree would contain the same leaves. This is because the method to select a valid cutting point are deterministic if only one dimension is available.
	 To avoid this, select randomly one element among available samples until it satisfies the validity conditions and use it as valid cutting point.
	*/
	if(nullptr == node.parent && 1==node.dim_cut.n_rows && valid_cut){

		dim=0;
		samples_dim = samples.col(dim);
		const double l = node.area(0,dim), r = node.area(1,dim);
		std::uniform_int_distribution<> dis(0,nMCsamples-1); // randomly select one row of current samples
		valid_cut = false;
		while(!valid_cut){
			cut=samples_dim(dis(gen)); // cutting point is randomly selected from samples along only feasible dimension				
			valid_cut = checkCut(cut,samples_dim,l,r); // cut need to be valid
		}
		// cut randomly picked but satisfies minimum length and minimum number of points
	}
	
};

// Method to establish if candidate cutting point cut is valid
bool Tree::checkCut(const double cut, const vec & samples_dim, const double l, const double r){
	
	bool valid_cut = false; // = true if and only if cut is a valid cutting point
	
	// check area condition for each block induced by cutting point cut
	if( ((cut-l)>= min_cut_length) && ((r-cut)>= min_cut_length) ){// not too small blocks
		
		valid_cut=true; // cutting point passed the condition on area
		
		// check number of MCMC samples along dimension dim that fall into each block of new partition
		const uvec left_idx=find(samples_dim<cut); // = {k in 0,...,samples_dim.n_rows-1 | samples_dim(k) < cut};
		const uvec right_idx=find(samples_dim>cut); // = {k in 0,...,samples_dim.n_rows-1 | samples_dim(k) > cut};
		
		// cut not valid if number of MCMC samples along dimension dim in each block is too small
		if( !( (left_idx.n_rows >= min_num_points_block) && (right_idx.n_rows >= min_num_points_block) ) ){
			valid_cut=false; 
		}
		
	}
	// valid_cut==false if cut did not satisfy AT LEAST ONE of the above conditions.
	
	return valid_cut;
};

// Method to determine the cutting point along a fixed dimension using Maximum Likelihood criterion
double Tree::ml_cutting_rule(const RawMCMC & rawMCMC, const Node & node, const uword dim){
	
	const uword M = rawMCMC.get_n_subsets();
	assert(M>0 && "Number of subsets must be at least 1!!!");
	
	const vec & samples_dim = rawMCMC.get_samples().col(dim);
	const vec & x = samples_dim(node.indices); // vector of samples in current node along dimension dim
	
	const uvec & mark = rawMCMC.get_mark()(node.indices);
	assert(all(mark>=0) && all(mark<M) && "!!! Indices of machines must be between 0 and M-1 !!!");
	
	const uvec idSorted_samples = sort_index(x); // index of sorted elements of x
	
	// Remove highest element in x: elements of partition have the following form [left,cut) [cut,right) hence do not consider highest value of current samples as candidate cut
	const uword num_cuts = x.n_rows-1; // number of candidate cutting points
	const uvec id_notMax = idSorted_samples.head(num_cuts); // ordered ID of candidate cutting points
	const vec cuts = x(id_notMax); // candidate cutting points: ordered values of samples, highest excluded
	const uvec sorted_mark = mark(id_notMax); // indices of subsets corresponding to ordered values of samples
	
	/* We need to determine the solution of the following optimization problem:
	 max{ \sum_{m=0}^{M-1} ( nLeft(m)*( \log(nLeft(m)) - \log(cut-l)) ) + \sum_{m=0}^{M-1} ( nRight(m)*( \log(nRight(m)) - \log(r-cut)) ) }
	 s.t.
	 		nLeft(m) + nRight(m) = nSubset(m)
	 		cut s.t. [l,cut) \cup [cut,r) = (l,r]
	 where
	 		nLeft(m) = number of samples in m-th subset falling into left area block		for m=0:M-1
	 		nRight(m) = number of samples in m-th subset falling into right area block		for m=0:M-1
	 		nSubset(m) = number of samples in m-th subset									for m=0:M-1	
	
	 For each candidate cutting point cuts(j) we determine the value of the objective function and then select the point j for which the value is maximum.	
	*/
	vec areaLeft = cuts - node.area(0,dim)*ones<vec>(num_cuts); // area of left block; left block corresponding to j-th cut: [l,cuts(j))
	uvec id_small = find(areaLeft<1e-4);
	areaLeft(id_small) = 1e-4*ones<vec>(id_small.n_rows); // avoid zeros
	
	vec areaRight = node.area(1,dim)*ones<vec>(num_cuts) - cuts; // area of right block; right block corresponding to j-th cut: [cuts(j),r)
	id_small = find(areaRight<1e-4);
	areaRight(id_small) = 1e-4*ones<vec>(id_small.n_rows);
	
	vec goodness_of_cuts; // objective value for each candidate cutting point
	goodness_of_cuts.zeros(num_cuts);
	for(uword m=0; m<M; m++){
		
		const uvec tmp_sub = find(mark==m);
		const uword nSubset = tmp_sub.n_rows; // number of samples (maximum included) belonging to m-th subset
		
		vec tmp;
		tmp.zeros(num_cuts);
		// if nSubset==0, first and second addend of objective function are 0 --> there's nothing to add to goodness_of_cuts!
		if(nSubset>0){
			const uvec cum_sum = cumsum(sorted_mark==m);
			for(uword j=0; j<num_cuts; j++){
				
				if(cum_sum(j)>0){ // consider first addend of objective function
					tmp(j) += cum_sum(j)*(std::log(cum_sum(j)) - std::log(areaLeft(j)));
				}
				
				if(cum_sum(j)<nSubset){ // consider second addend of objective function
					tmp(j) += (nSubset - cum_sum(j))*(std::log(nSubset - cum_sum(j)) - std::log(areaRight(j)));
				}
				
			}				
		}
		
		goodness_of_cuts += tmp; // update objective value summing result of m-th subset
	}

	const uword index = id_notMax(index_max(goodness_of_cuts)); // index of sample corresponding to ML cutting point	
	
	return x(index);

};

// Method to grow the binary tree
void Tree::grow(const RawMCMC & rawMCMC, const uvec & dimCut, const bool verbose, std::vector<MCestimate> & leaves){
	
	Node node; // create root node
	node.area = rawMCMC.get_area();
	const uword nMCsamples = rawMCMC.get_samples().n_rows; // number of rows in samples matrix
	assert(nMCsamples>0 && "Number of samples must be at least 1!!");
	node.indices = linspace<uvec>(0,nMCsamples-1,nMCsamples); // at first step consider all rows of RawMCMC object
	node.dim_cut = dimCut;
	
	build(rawMCMC, node, verbose, leaves); // grow binary tree
	
	// additional normalization of tree ONLY if more than 1 subset involved
	if(rawMCMC.get_n_subsets()>1) normalize(leaves);
	
};

// Method to recursively grown the binary tree
void Tree::build(const RawMCMC & rawMCMC, Node & node, const bool verbose, std::vector<MCestimate> & leaves){
	
	const uword nMCsamples = node.indices.n_rows; // number of samples in current node
	assert(nMCsamples>0 && "Number of samples must be at least 1!!");
	
	const mat & area = node.area; // select area corresponding to current node
	if(verbose){
		area.print("Sampling area at current node: ");
	}
	
	bool valid_cut = false; // true if valid partition found
	
	// Initialize objects for cut search
	double cut=0.0; // cutting point
	uword dim=0; // dimension along which cut is searched
	vec samples_dim; // samples' matrix reduced to column dim

	// OBSERVE: Only at root node we could have nMCsamples<=min_num_points_block since it could happen that rawMCMC stores less elements than min_num_points_block.
	if(nullptr!=node.parent){
		assert(nMCsamples>=min_num_points_block && "Number of samples must be at least min_num_points_block at non-root nodes!!");
	}
	
	// if current node stores less samples than min_num_points_block and has no candidate dimensions along which cutting point should be searched, the node is automatically leaf!
	if(nMCsamples>min_num_points_block && node.dim_cut.n_rows>0){
		findCut(rawMCMC,node,cut,dim,samples_dim,valid_cut); // search for valid cut
	}
	
	if(valid_cut){ // Valid partition has been found --> define the partition to be associated to children nodes
		
		assert(nMCsamples>=min_num_points_block && "Number of samples in node is less than minimum!!!");
		
		if(verbose) std::cout << "\n Valid partition found! Cutting point " << cut << " along dimension " << dim << " is valid! Continue searching partition into children nodes..." << std::endl << std::endl;		
		
		/* Randomly split samples that coincide with cutting point in left and right area (helps when sampling method get stucked) */
		// Random devices
		std::random_device rd;
		std::mt19937 gen(rd());
		std::bernoulli_distribution dis(0.5);
		uvec left_idx, right_idx;
		if(dis(gen)){ // add points coinciding with cutting point to right partition's block
			right_idx=find(samples_dim>=cut);
			left_idx=find(samples_dim<cut);
		}else{ // add points coinciding with cutting point to left partition's block
			right_idx=find(samples_dim>cut);
			left_idx=find(samples_dim<=cut);
		}		
		assert( left_idx.n_rows+right_idx.n_rows==nMCsamples && "!!!Something wrong in splitting samples into left and right!!!");		
		assert( left_idx.n_rows>=min_num_points_block && "Number of samples falling left w.r.t. valid cutting point is lower than minimum!!");
		assert( right_idx.n_rows>=min_num_points_block && "Number of samples falling right w.r.t. valid cutting point is lower than minimum!!");
		
		// create children of current node
		Node nodeLeft, nodeRight; // created as root nodes (nullptr for parent, childLeft, childRight)
		node.generateChild(dim,cut,min_cut_length,true,left_idx,nodeLeft); // defines left child
		node.generateChild(dim,cut,min_cut_length,false,right_idx,nodeRight); // defines right child
		
		build(rawMCMC, nodeLeft, verbose, leaves); // proceed search of valid partition on left child
		build(rawMCMC, nodeRight, verbose, leaves); // proceed search of valid partition on right child
		
	}else{ // No valid partition has been found --> this is one of the element of the partition!!

		if(verbose)	std::cout << "\n Leaf node reached!" << std::endl << std::endl;
		
		leaves.emplace_back(rawMCMC,node.indices,area); // store into vector of leaves statistics on samples falling in current node 
		
		if(verbose){
			leaves.back().mean.print("Mean");
			leaves.back().Cov.print("Cov");
			leaves.back().sampling_area.print("Area");
			std::cout << "Log_prob " << leaves.back().log_prob << std::endl;
		}
		
	}
	
};

// Method to normalize the binary tree (represented as leaf nodes)
void Tree::normalize(std::vector<MCestimate> & leaves){
	
	const unsigned int n_leaves = leaves.size();
	
	// Logarithm of each prob_hat(k) as computed at each leaf node:
	// \log(prob_hat(k)) := -(M-1) * \log(|block(k)|) + \sum_{m=0}^{M-1} ( \log(n_block(m,k)) - \log(N(m)) )
	vec logProbHat(n_leaves);
	for(unsigned int i=0; i<n_leaves; i++){
		logProbHat(i) = (leaves[i]).log_prob;
	}
	
	// Logarithm of normalizing constant: \log(Z) = \log( \sum_{k=0}^{K-1} \exp( \log(prob_hat(k)) ) )
	// (to avoid too high values, substract the maximum)
	double max_log_prob = max(logProbHat);	
	double log_normalizer = max_log_prob + std::log(accu(exp(logProbHat-max_log_prob*ones<vec>(n_leaves))));
	
	// Store prob(k) values in logarithmic form: \log(prob(k)) = \log(prob_hat(k)) - \log(Z)
	for(auto & i: leaves){
		i.log_prob -= log_normalizer;
	}

};