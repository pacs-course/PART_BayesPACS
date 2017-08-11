//#define ARMA_NO_DEBUG // if defined, avoid all checks of dimensions through ()
//#define ARMA_DONT_PRINT_ERRORS // if defined, errors and warnings not print on screen
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <chrono>
#include <string>
#include <armadillo>
#include <vector>
#include "Parameters.hpp"
#include "aggregatedResampling.hpp"

using namespace arma;

int main(int argc, char** argv){
	
	const auto t_start = std::chrono::high_resolution_clock::now();
	
	if(argc<3 || argc>4){ // At least input and output paths, at most input, output and parameter paths
		std::cerr << "!!! ERROR: Program requires two or three inputs !!!" << std::endl;
		std::cerr << "First input: paths for input subchains. " << std::endl; 
		std::cerr << "Second input: path for PART output. " << std::endl;
		std::cerr << "Third input (OPTIONAL): path of parameters' file. " << std::endl;
		return 1;
	}
	
	// at least 2 inputs, at most 3 inputs
	const std::string pathInputDataFile = argv[1], outputFile=argv[2];
	std::string paramsFile;
	if(argc==3){ // Use default parameters
		std::cout << "Use default parameters defined in parameters/defaultParams.par file."<< std::endl;
	}else{
		paramsFile = argv[3];
		std::cout << "Path of parameters' file: " << paramsFile << std::endl;		
	}

	try{

		// read parameters for PART from file
		Parameters par(paramsFile); // Initializes N, d, n_samples to zeros
		
		// read path of each subchain file
		std::ifstream fileInput (pathInputDataFile);
		if(!fileInput.is_open()){
			throw std::logic_error ("Unable to open file of paths!!");
		}
		
		std::vector<std::string> pathSubchain(par.M);
		uword pos=0;
		while(!fileInput.eof()){
			std::getline(fileInput,pathSubchain[pos]);
			pos++;
		}
		fileInput.close();
		
		if(par.M!=pathSubchain.size()){
			throw std::length_error ("Different number of subsets specified in parameters file!!");
		}

		// store subchains
		std::vector<mat> subchains;
		subchains.reserve(par.M);
		
		uword dim_data = 0;
		for(unsigned int m=0; m<par.M; m++){ // select groups of data
			
			mat data;
			if(!data.load(pathSubchain[m])){ //
				throw std::logic_error ("Unable to open data file!!");
			}
			
			const uword n_data = data.n_rows;
			const uword dim_inner = data.n_cols;
			
			if(dim_inner<1 || n_data<1){
				throw std::length_error ("Each input matrix must contain at least 1 row and 1 column !!!");
			}

			if(m==0){
				dim_data=dim_inner;
			}

			if(dim_inner != dim_data){
				throw std::logic_error ("Number of parameters to be monitored do not coincide with parameters' one!!");
			}
			
			subchains.push_back(data);
			
			par.N(m) = n_data;
			
			std::cout << "Dimension of subchain: Number of samples: " << n_data << " Number of features: " << dim_data << std::endl;
			
		}
		
		par.d = dim_data; // Number of features
		
		par.n_samples = std::max(max(par.N),par.resample_N); // Number of MCMC samples to be drawn from the aggregated posterior
		
		par.display(); // display basic information on parameters
		
		// perform the posterior aggregation on each subset and draw samples from the estimated posterior 
		mat aggr_post_samples;
		aggregatedResampling(par, subchains, aggr_post_samples); 
		
		// save on file samples drawn from aggregated posterior
		if(!aggr_post_samples.save(outputFile, csv_ascii)){
			throw std::logic_error ("Cannot save results in output file!!");
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