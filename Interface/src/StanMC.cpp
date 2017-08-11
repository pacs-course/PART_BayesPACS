#include <iostream>
#include <chrono>
#include <stdexcept>
#include <cassert>
#include <string>
#include <vector>
#include "StanMC.hpp"

using namespace arma;

// Constructor from file
StanParams::StanParams(const std::string & fileMC){

	// read path of each subchain file
	std::ifstream fileInput (fileMC);
	if(!fileInput.is_open()){
		throw std::logic_error ("Unable to open file of paths!!");
	}
	
	std::vector<std::string> inputMC;
	std::string tmpStr;
	while(std::getline(fileInput,tmpStr)){
		inputMC.push_back(tmpStr);
	}
	fileInput.close();
	
	assert(6==inputMC.size() && "Wrong number of input parameters!!");
	
	// std::stol : convert to long int
	// std::stoi : convert to int
	m_samplingIter = static_cast<uword>(std::stol(inputMC[0]));
	m_burnin = static_cast<uword>(std::stol(inputMC[1]));
	m_thin = static_cast<uword>(std::stol(inputMC[2]));
	
	m_pathR = inputMC[3];
	m_pathExeStan = inputMC[4];
	
	m_useInit = static_cast<bool>(std::stoi(inputMC[5]));
	
};

// Split data and generate M independent suchains
void StanMC::createInputSubchains(const std::string & pathInputData, const uword M, MCMCParams * parMC){

	// dynamic cast of MCMCParams
	StanParams * ptrParStan = dynamic_cast<StanParams*>(parMC);
	assert(nullptr!=ptrParStan && "Dynamic cast failed!!");
	
	// store parameters as attribute
	m_parStan = *ptrParStan;
	
	// load input dataset
	mat data;
	if(!data.load(pathInputData)){
		throw std::logic_error ("Unable to read data file!!");
	}
	const uword dim = data.n_cols;
	const uword n_data = data.n_rows;
	if(dim<1 || n_data<1){
		throw std::length_error ("Input matrix must contain at least 1 row and 1 column !!");
	}
	
	std::cout << "Dimension of dataset: Number of samples = " << n_data << "; Number of features = " << dim << std::endl;
	
	if(n_data<=M){
		throw std::logic_error ("Impossible split: number of subsets to be created larger than number of rows of input matrix !!");
	}
	
	// number of rows in each subdata (except last subdata that will contain all remaining rows in data matrix)
	const uword nRowsSubset = static_cast<uword>(std::floor(n_data/static_cast<double>(M)));
	if(nRowsSubset<1){
		throw std::length_error ("Subdata matrix cannot contain less than one row!!");
	}
	
	// Reserve size for container of chains
	m_subchains.reserve(M);

	// split dataset: all subsets have nRowsSubset rows, except last one that stores all remaining rows
	for(unsigned int m=0; m<M-1; m++){
		
		const auto t_start = std::chrono::high_resolution_clock::now();
		
		// store m-th subdata
		mat subdata(data.rows(0+m*nRowsSubset,(m+1)*nRowsSubset-1));
		
		// generate m-th subchain and store value into m_subchains[m]
		generateSubchain(subdata,m);
		
		const auto t_end = std::chrono::high_resolution_clock::now();
		
		std::cout << "Execution time to generate " << m << "th chain : " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count() << " milliseconds." << std::endl;
		
	}
	
	const auto t_start = std::chrono::high_resolution_clock::now();	

	// save all remaining rows
	mat subdata(data.rows((M-1)*nRowsSubset,data.n_rows-1));
		
	// generate last subchain and store value into m_subchains[M-1]
	generateSubchain(subdata, M-1);

	const auto t_end = std::chrono::high_resolution_clock::now();
	
	std::cout << "Execution time to generate " << M << "th chain : " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count() << " milliseconds." << std::endl;
	
};

// Generate one subchain indexChain using CmdStan
void StanMC::generateSubchain(const mat & subdata, const unsigned int indexChain){
	
	// path of R files
	const std::string pathR = m_parStan.get_pathR();
	
	// save subdata onto .csv file
	const std::string nameSubdata(pathR + "/subdata.csv");
	if(!subdata.save(nameSubdata, csv_ascii)){
		throw std::logic_error ("!!! ERROR saving subdata !!! ");
	}
	
	// transform data into .data.R extension using R
	const std::string cmdR_string ("Rscript "+ pathR + "/dataDump.R");
	const char * cmdR = cmdR_string.c_str();
	if(std::system(cmdR)!=0){
		throw std::logic_error ("!!! ERROR executing R session!!! ");
	}
	// content of pathR+"/subdata.csv" stored in .data.R extension
	
	// Run CmdStan 
	std::string cmdStan_string(m_parStan.get_pathExeStan() + " sample num_samples=" + std::to_string(m_parStan.get_samplingIter())
								+ " num_warmup=" + std::to_string(m_parStan.get_burnin()) + " thin=" + std::to_string(m_parStan.get_thin())
								+ " data file=" + pathR + "/subdata.data.R");
	
	// if user wants to initialize variables for MC, generate .init.R and add it to CmdStan command
	if(m_parStan.get_useInit()){
		const std::string cmdRInit_string ("Rscript "+ pathR + "/initDump.R");
		const char * cmdRInit = cmdRInit_string.c_str();
		if(std::system(cmdRInit)!=0){
			throw std::logic_error ("!!! ERROR executing R session!!! ");
		}
		cmdStan_string = cmdStan_string + " init=" + pathR + "/inits.init.R";		
	}
	
	// path default output file of CmdStan
	const std::string nameDefaultChain (pathR + "/output.csv");
	
	cmdStan_string = cmdStan_string + " output " + nameDefaultChain;
	const char * cmdStan = cmdStan_string.c_str();
	if(std::system(cmdStan)!=0){
		throw std::logic_error ("!!! ERROR executing CmdStan session!!! ");
	}
	// generate subchain in pathR+"/output.csv"
		
	// Remove comments (#) and header line and store result into .csv file
	const std::string nameSubchain (pathR + "/subchain" + std::to_string(indexChain) + ".csv");
	const std::string cmdGrep_string("grep \"#\" -v "+ nameDefaultChain + " | grep -v lp__ > " + nameSubchain);
	const char * cmdGrep = cmdGrep_string.c_str();
	if(std::system(cmdGrep)!=0){
		throw std::logic_error ("!!! ERROR removing comments and header from subchain file!!! ");
	}

	// Remove subdata.csv, .data.R, .init.R and temporary output files
	const std::string cmdRemove_string("rm " + nameSubdata + " " + pathR + "/subdata.data.R " + pathR + "/inits.init.R " + nameDefaultChain);
	const char * cmdRemove = cmdRemove_string.c_str();
	if(std::system(cmdRemove)!=0){
		throw std::logic_error ("!!! ERROR removing subdata, inits or temporary subchain file!!! ");
	}
	
	mat chain;
	if(!chain.load(nameSubchain)){
		throw std::logic_error ("!!! ERROR loading subchain !!! ");
	}
	
	m_subchains.push_back(chain);
	
};