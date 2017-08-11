#include <iostream>
#include <stdexcept>
#include "Parameters.hpp"
#include "sigpack/sigpack.h" // read input parameters

Parameters::Parameters(const std::string & paramsFile){
	
	sp::parser testpar(paramsFile); // Create parser to read from input file (already checks if not able to open file!)
	
	M = testpar.getParam<uword>("M", 4); // Number of subsets to be used
	if(M<1){
		throw std::length_error ("Too few input subsets!!");
	}
	
	pairwise_aggregation = testpar.getParam<bool>("pairwise_aggregation", 0); // Boolean of pairwise aggegation
	improve_matching = testpar.getParam<bool>("improve_matching", 1); // Boolean of Improve Matching
	halve = testpar.getParam<bool>("halve", 0); // Boolean of Halve
	
	ntree = testpar.getParam<uword>("ntree", 16); // Number of trees
	if(ntree<1){
		throw std::length_error ("At least 1 tree required to build the forest!!");	
	}
	
	kd_cut = testpar.getParam<bool>("kd_cut", 1); // Boolean of KD-Tree Partition 
    
	min_frac_points_block = testpar.getParam<double>("min_frac_points_block", 0.01); // Minimum fraction of points per block
	
	min_cut_length = testpar.getParam<double>("min_cut_length", 0.001);
	
	resample_N = testpar.getParam<uword>("resample_N", 10000); // Number of sample to be resampled
	if(resample_N<1){
		throw std::length_error ("Too few number of samples to be drawn from the estimated posterior!!");	
	}
	
	gaussian_smooth = testpar.getParam<bool>("gaussian_smooth", 1); // Boolean of Gaussian Smoothing
	
	verbose = testpar.getParam<bool>("verbose", 0); // Verbose
	
	N.zeros(M);
	
	d = 0;
	
};

void Parameters::display() const {
	
	std::cout << "*************** Information on input parameters ***************" << std::endl;
    
	std::cout << "Number of features: " << d << std::endl;
	
	std::cout << "Number of subsets: " << M << std::endl;
	N.print("Number of MCMC samples stored on each subset:");
	
	if(pairwise_aggregation){
		std::cout << "Use pairwise aggregation with";
		if(improve_matching){
			std::cout << " improved matching for subsets' indices combination ";
		}else{
			std::cout << " simple matching for subsets' indices combination ";
		}
		if(halve){
			std::cout << "halving the minimum fraction of MCMC samples in partition's blocks at each aggregation stage." << std::endl;
		}else{
			std::cout << "without halving the minimum fraction of MCMC samples in partition's blocks at each aggregation stage." << std::endl;
		}
	}else{
		std::cout << "Use one-stage aggregation. Ignore improve_matching and halve parameters ..." << std::endl;
	}
	
	std::cout << "Build a forest of "<< ntree << " trees. " << std::endl;
	
	if(kd_cut){
		std::cout << "Use median/KD-method to determine cutting points." << std::endl;
	}else{
		std::cout << "Use ML-method to determine cutting points." << std::endl;
	}
	
	if(min_frac_points_block>0.1){
		std::cout << "WARNING: Minimum fraction of samples in each block " << min_frac_points_block << " is quite high. Consider reducing to improve estimates ! " << std::endl;
	}else{
		std::cout << "Minimum fraction of MCMC samples in each partition's block: " << min_frac_points_block << std::endl;
	}
	
	std::cout << "Minimum length of each partition's block along each dimension:" << min_cut_length << std::endl;
	
	std::cout << "Effective number of MCMC samples to be drawn from the aggregated posterior: n_samples = max(max(N),resample_N) = " << n_samples << std::endl;
	
	if(gaussian_smooth){
		std::cout << "Use local Gaussians to approximate the aggregated posterior. " << std::endl;
	}else{
		std::cout << "Use local uniforms to approximate the aggregated posterior. " << std::endl;
	}

	if(verbose){
		std::cout << "All messages in code will be displayed." << std::endl;
	}else{
		std::cout << "Only relevant information will be displayed."<< std::endl;
	}
	
};