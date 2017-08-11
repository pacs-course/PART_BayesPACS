#ifndef H_MCMC_HPP
#define H_MCMC_HPP

#include <string>
#include <vector>
#include <armadillo>
using namespace arma;

//! Base struct for parameters of MCMC method
struct MCMCParams{
	virtual void f(){};
};

//! Base class for MCMC generation
class MCMC{
	
	protected:
		
		std::vector<mat> m_subchains; /*!< vector of matrices; each matrix stores MCMC samples from one subset. */

		//! virtual method; generate subchain with index indexChain based on subdata.
		virtual void generateSubchain(const mat & subdata, const unsigned int indexChain) = 0;
		
	public:
		
		//! virtual method; generate all subchains.
		virtual void createInputSubchains(const std::string & pathInputData, const uword M, MCMCParams * parMC) = 0;
		
		//! get method for subchains attribute
		inline const std::vector<mat> & get_subchains() const { return m_subchains; }
		
};
#endif