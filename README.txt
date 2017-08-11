README for C++11 implementation of PART algorithm.

Reference paper (Wang et alt.) and MATLAB implementation available at: https://github.com/richardkwo/random-tree-parallel-MCMC

----------------------------------------------------------------------------------------------------------------------------------------

PRELIMINARY INFORMATION:

1) Install ARMADILLO library http://arma.sourceforge.net/download.html
	OBSERVE: we used version 7.800.2 of Armadillo library.

2) (OPTIONAL) Install following programs for:
	- documentation: Doxygen http://www.stack.nl/~dimitri/doxygen/
	- profiling: Gprof https://sourceware.org/binutils/docs/gprof/
	- memory usage: Valgrind http://valgrind.org/

3) Set environmental variables in setEnv.sh file and the source:
	> source setEnv.sh

4) To read PART parameters from file, we use SigPack. It is not required to download SigPack since we already include it (see ./PART/include/sigpack).
   For additional information on how to read parameters from file using SigPack see http://sigpack.sourceforge.net/

5) Beside PART algorithm, we also developed an interface (./Interface) to generate MCMC samples from subsets of data.
   The user could use the desired generator, defining accordingly a derived class from ./Interface/include/MCMC.hpp.
   As generator we choose CmdStan, that is already available in ./Interface/cmdstan-2.16.0. For additional details on usage, please refer to http://mc-stan.org/users/interfaces/cmdstan.
   To make use of our implementation of PART interface, the user should also install R (https://cran.r-project.org/) and RStan (https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux), used to transform data in compatible format with respect to CmdStan.
   
   OBSERVE: CmdStan execution may require long time, according to dimensions of input dataset and to the model complexity. Because of this we keep separate folders for Interface and PART, to easily run synthetic tests on already generated subchains.
----------------------------------------------------------------------------------------------------------------------------------------

COMPILE AND RUN - PART
Instructions for PART implementation only (./PART).

1) Makefile options:
	
	To compile our program in optimization mode (using -O2 flag):
	> make

	To compile our program in debug mode (using -g -O0 flags):
	> make debug

	To compile our program in profile mode (-pg flag):
	> make prof
	
	To generate Doxygen documentation:
	> make doxy

	To remove objects and executable:
	> make clean

	To restore folder:
	> make distclean


2) Execution

	Our program requires at least two input arguments to run, an could use an optional additional input:
    1. Path of .txt file into which the paths of files storing subchains have been defined. Each subchain's file should have a compatible extension with respect to \textit{Armadillo} (see \url{http://arma.sourceforge.net/docs.html#save_load_mat}).
    2. Path of output file (.csv, no headers; see csv_ascii at http://arma.sourceforge.net/docs.html#save_load_mat), storing MCMC samples drawn from the aggregated posterior obtained through PART algorithm.
	3. (OPTIONAL) Path of PART parameters file, defining parameters to be used in PART algorithm. This is a .par file for compatibility with SigPack (see http://sigpack.sourceforge.net/ for additional details).
		In case this file is not given as input to the program, default parameters will be used, as defined ./examples/parameters/PART/defaultParamsPART.par.
		Recall that all parameters defined in this file must be coherent with the input subchains in order for the program to properly run.

----------------------------------------------------------------------------------------------------------------------------------------
Examples - PART:

Compile and Run - simple:
	> cd PART
	> make
	> ./main ../examples/data/d1/chains/M10/input_M10_N5.txt ../examples/data/d1/outPART_M10_N5_kdPairSmooth.csv ../examples/parameters/PART/M10/kdPairSmooth.par

Compile and Run - synthetic tests (see Report_Raciti_Riva.pdf)
To run all synthetic test for a fixed dimension d, a bash script is already available in each data folder.
E.g. from PART directory , to run all examples for dimension d=9, proceed as follows:
	> chmod 777 ../examples/data/d9/runPART.sh
	> ../examples/data/d9/runPART.sh

----------------------------------------------------------------------------------------------------------------------------------------

COMPILE AND RUN - PART + Interface with CmdStan
Instructions for MCMC generation of suchains using CmdStan, to be used as input of PART.
Please refer to ./Interface directory.

1) Makefile options:
	
	To compile our program in optimization mode (using -O2 flag):
	> make

	To compile our program in debug mode (using -g -O0 flags):
	> make debug

	To compile our program in profile mode (-pg flag):
	> make prof
	
	To generate Doxygen documentation:
	> make doxy

	To remove objects and executable:
	> make clean

	To restore folder:
	> make distclean

2) Execution

	PART program interfaced with CmdStan, requires at lest three inputs to run, an could use an optional additional input:
	1. Path of data file into which the entire dataset is stored. To read input data we use load method of Armadillo Mat class, that should automatically recognize file extensions (see http://arma.sourceforge.net/docs.html#save_load_mat).
	2. Path of PART output, storing MCMC samples from aggregated posterior in .csv extension (see csv_ascii at http://arma.sourceforge.net/docs.html#save_load_mat).
	3. Path of .txt file storing information on MCMC parameters.
	   In our case, this file must contain all information to run CmdStan and R, hence the following ordered inputs:
			samplingIter: number of sampling iterations of MCMC algorithm;
			burnin: number of burnin (warmup) iterations of MCMC algorithm;
			thin: thinning paramter for MCMC algorithm;
			pathR: path into which R scripts are located (i.e. working directory defined inside R scripts);
			pathExeStan: path into which executable of Stan model is located (i.e. the one obtained compiling .stan model using CmdStan; see Report_Raciti_Riva.pdf);
			useInit: flag to establish if input file for initializations of parameters for MCMC sampling should be used; if set to true, the file pathR/inits.init.R will be used.
	   An example of the described file storing default values for CmdStan interface is available at ./examples/parameters/MCMC/defaultParamsMCMC.txt
	4. (OPTIONAL) Path of PART parameters file, defining parameters to be used in PART algorithm. This is a .par file for compatibility with SigPack (see http://sigpack.sourceforge.net/ for additional details).
		In case this file is not given as input to the program, default parameters will be used, as defined ./examples/parameters/PART/defaultParamsPART.par.
		Recall that all parameters defined in this file must be coherent with the input subchains in order for the program to properly run.

----------------------------------------------------------------------------------------------------------------------------------------
Examples - PART + Interface:

Compile and Run - simple:
To compile and run a simple example, proceed as follows:
	> cd Interface
	> make clean
	> make
	> ./main ../examples/data/d1/data.csv ../examples/data/d1/outPART_M10_N5_kdPairSmooth.csv ../examples/parameters/MCMC/paramsMCMC_N5.txt ../examples/parameters/PART/M10/kdPairSmooth.par

Compile and Run - all PART configurations:
To compile and run all configurations for a fixed dimension d, run the bash script located in the associated folder.
E.g. for dimension d=9, proceed as follows:
	> cd Interface
	> make clean
	> make
	> chmod 777 ../examples/data/d9/runInterface.sh
	> ../examples/data/d9/runInterface.sh

----------------------------------------------------------------------------------------------------------------------------------------

Observe: to generate Doxygen documentation in doc directory, move to PART directory or to Interface directory and type:
> make doxy

----------------------------------------------------------------------------------------------------------------------------------------

CODE STRUCTURE

PART_BayesPacs_master 									Master directory
|		
|		
|---> README.txt										README file								 
|
|
|---> setEnv.sh											shell file to set environmental variables
|
|
|---> doc												Directory for documentation
|		|
|		|---> Doxyfile										Doxygen file to generate documentation.
|		|
|		|---> PART-Wang15.pdf								Original paper
|		|
|		|---> Report_Raciti_Riva.pdf						Project's report
|
|
|
|---> PART												Directory of C++11 implementation of PART algorithm
|		|
|		|---> Makefile										Makefile
|		|
|		|---> src											Source code for PART only
|		|
|		|---> include										Headers for PART only
|				|
|				|---> sigpack									SigPack headers
|
|
|
|---> Interface											Directory for interface code
|		|
|		|---> Makefile										Makefile
|		|
|		|---> src											Source code for Interface
|		|
|		|---> include										Headers for Interface
|		|
|		|---> cmdstan-2.16.0								CmdStan directory (already built)
|		|
|		|---> Rdir											Directory of R scripts
|
|
|
|---> examples											Directory for synthetic tests
		|
		|---> parameters
		|		|
		|		|---> MCMC									Directory for MCMC parameters' files
		|		|
		|		|---> PART									Directory for PART parameters' files
		|
		|
		|---> data										Directory storing data for different dimensions d
		|		|
		|		|---> d1
		|		|
		|		|---> d9
		|		|
		|		|---> d29
		|
		|---> MATLAB code								Directory for MATLAB files (accuracy)
