#ifndef H_STANMC_HPP
#define H_STANMC_HPP

#include <string>
#include <armadillo>
#include "MCMC.hpp"

using namespace arma;

//! Struct for parameters of MCMC method using Stan.
struct StanParams : public MCMCParams{
	
	uword m_samplingIter; /*!< Number of sampling iterations of MCMC method. */
	uword m_burnin; /*!< Number of burnin iterations of MCMC method. */
	uword m_thin; /*!< Thinning of MCMC method. */
	
	std::string m_pathR; /*!< Path into which all R files are located. */
	std::string m_pathExeStan; /*!< Path into which executable of Stan model is located */
	
	bool m_useInit; /*!< \c true if user decides to define initialization file for MCMC generation. */
	
	public:
		
		//! Empty constructor
		StanParams(){};
		
		//! Constructor from file
		StanParams(const std::string & fileMC);
		
		//! get method for m_samplingIter
		inline const uword get_samplingIter() const { return m_samplingIter; }
		
		//! get method for m_burnin
		inline const uword get_burnin() const { return m_burnin; }
		
		//! get method for m_thin
		inline const uword get_thin() const { return m_thin; }		
		
		//! get method for m_pathR		
		inline const std::string & get_pathR() const { return m_pathR; }
		
		//! get method for m_pathExeStan
		inline const std::string & get_pathExeStan() const { return m_pathExeStan; }
		
		//! get method for m_useInit
		inline const bool get_useInit() const { return m_useInit; }
		
};

//! Class for MCMC generation through Stan.
class StanMC: public MCMC{
	
	StanParams m_parStan; /*!< parameters for MCMC sampling using Stan. */
	
	//! Generate one subchain indexChain using CmdStan
	/*!
	@param subdata (\c mat): subset of data to be used to generate chain.
	@param indexChain (\c unsigned int): index of subchain currently generated.
	*/
	void generateSubchain(const mat & subdata, const unsigned int indexChain);
	
	public:
		
		//! Empty constructor
		StanMC(){};
		
		//! Split data and generate M independent suchains
		/*!
		@param pathInputData (\c string): path of entire dataset.
		@param M (\c uword): number of subchains to be generated.
		@param parMC (\c MCMCParams \c * ): pointer to MCMCParams object, storing information on how to generate chains.
		
		@par Given the entire dataset, this method splits it into M subdata and uses each of them to generate M subchains independently.
		The result will be stored in private attribute m_subchains.
		
		@see MCMC basic class.
		@see generateSubchain for details on how to generate one single subchain using Stan.
		*/
		void createInputSubchains(const std::string & pathInputData, const uword M, MCMCParams * parMC);
		
		//! get method for parameters' attribute.
		inline const StanParams & get_parStan() const { return m_parStan; }
		
};
#endif