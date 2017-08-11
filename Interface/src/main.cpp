//#define ARMA_NO_DEBUG // if defined, avoid all checks of dimensions through ()
//#define ARMA_DONT_PRINT_ERRORS // if defined, errors and warnings not print on screen
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <chrono>
#include <string>
#include <memory>
#include <new>
#include <armadillo>
#include <vector>
#include "Parameters.hpp"
#include "StanMC.hpp"
#include "aggregatedResampling.hpp"

using namespace arma;

int main(int argc, char** argv){
	
	const auto t_start = std::chrono::high_resolution_clock::now();

	// check correctness of input
	if(argc<4 || argc>5){
		std::cerr << "!!! ERROR: Program requires three or four inputs !!!" << std::endl;
		return 1;
	}
	
	const std::string pathInputData = argv[1], outputAggrFile=argv[2], pathFileMC=argv[3];
	
	std::string paramsFile; // optional input
	if(argc==3){ // Use default parameters and no initialization!!!
		std::cout << "Use default parameters."<< std::endl;
	}else{
		paramsFile = argv[4];
	}
	
	try{

		// read parameters for PART algorithm from file
		Parameters par(paramsFile);
		
		// read parameters for MCMC method from file
		StanParams parStan(pathFileMC);
		
		// create unique pointer to derived class StanMC
		std::unique_ptr<StanMC> ptrMCStan (new StanMC);
		
		// call method to split data and generate subchains
		ptrMCStan->createInputSubchains(pathInputData,par.M,&parStan);
		
		// store subchains to be used as input of PART algorithm
		std::vector<mat> subchains(ptrMCStan->get_subchains());
		
		assert(par.M==subchains.size() && "ERROR: number of chains generated differ from desired one!!");
		
		// Set other attributes of Parameters object
		par.d = subchains.begin()->n_cols;
		
		const uword N_value = parStan.get_samplingIter()/static_cast<double>(parStan.get_thin());
		if(N_value<1 || par.d<1){
			throw std::length_error ("Each subchain matrix must contain at least 1 row and 1 column !!");
		}
		par.N = N_value * ones<uvec>(par.M); // Number of samples stored on each subchain
		
		par.n_samples = std::max(max(par.N),par.resample_N); // Number of MCMC samples to be drawn from the aggregated posterior
		if(par.n_samples<1){
			throw std::length_error ("Cannot resample less than 1 sample !!");
		}
		
		par.display(); // display basic information on parameters

		// perform the posterior aggregation on each subset and draw samples from the estimated posterior 
		mat aggr_post_samples;
		aggregatedResampling(par, subchains, aggr_post_samples); 
		
		// save on file samples drawn from aggregated posterior
		if(!aggr_post_samples.save(outputAggrFile, csv_ascii)){
			throw std::logic_error ("!!! ERROR saving results in output file !!! ");
		}
		
	}
	catch(const std::length_error& length_e){
		std::cerr << "Length error: " << length_e.what() << std::endl;
	}
	catch(const std::logic_error& logic_e){
		std::cerr << "Logic error: " << logic_e.what() << std::endl;
	} 
	
	const auto t_end = std::chrono::high_resolution_clock::now();
	std::cout << "Execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count() << " milliseconds." << std::endl;

	return 0;
	
};