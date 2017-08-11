#ifndef H_RAWMCMC_HPP
#define H_RAWMCMC_HPP

#include <armadillo>
#include <vector>

using namespace arma;

/*! @file RawMCMC.cpp
@brief Basic class to handle information on MCMC samples drawn on different subsets.
*/
class RawMCMC{

	public:

		//! Copy constructor
		RawMCMC(const RawMCMC & r)= default;

		//! Copy assignment operator
		RawMCMC& operator=(const RawMCMC & r)= default;
		
		//! Constructor from vector of pointers to matrices of MCMC samples drawn on each subset
		/*!
		@param subchains (\c vector \c < \c const \c mat \c * \c >): vector of pointers to matrices of MCMC samples.

		To create \c m_area (representing sampling area), call createDefaultArea() which computes sampling area from \c m_samples.

		To create \c m_samples and \c m_mark:\n
		for each subset m=0:subchains.size()-1\n
		- in \c m_samples concatenate MCMC samples corresponding to subset m,
		- in \c m_mark store index m of subset from which samples have been drawn.
		*/
		RawMCMC(const std::vector<const mat*> & subchains);

		//! getter of \c m_M attribute
		const uword get_n_subsets() const;

		//! getter of \c m_sumLogN attribute
		const double get_sumLogN() const;
		
		//! getter of \c m_samples attribute
		const mat & get_samples() const;
		
		//! getter of \c m_mark attribute
		const uvec & get_mark() const;
		
		//! getter of \c m_area attribute
		const mat & get_area() const;

		//! setter of \c m_M attribute
		void set_n_subsets(const uword nSubsets);

		//! setter of \c m_sumLogN attribute
		void set_sumLogN(const double sum_log_N);
		
		//! setter of \c m_samples attribute
		void set_samples(const mat & samples);
		
		//! setter of \c m_mark attribute
		void set_mark(const uvec & mark);
		
		//! setter of \c m_area attribute
		void set_area(const mat & area);
		
	private:
		
		static constexpr double enlargeAreaFactor = 1.001; /*!< Factor to enlarge sampling area. */
		
		uword m_M; /*!< Number of subsets. */
		double m_sumLogN; /*!< Sum of logarithm of number of rows for each subset's matrix. */
		mat m_samples; /*!< matrix of MCMC samples drawn from each subsets, all stored together. */
		uvec m_mark; /*!< Vector of subset's indices for each row in m_samples; used to retrieve from which subset the MCMC samples have been drawn. */
		mat m_area; /*!< Sampling area, represented as left and right extreme along each dimension. */
		
		//! Method to compute and set sampling area from MCMC samples
		/*!
		The sampling area is computed as follows:\n
		for each dimension p=0:m_samples.n_cols-1\n
			- area(0,p) represents the minimum boundary (left extreme) along dimension p\n
			- area(1,p) represents the maximum boundary (right extreme) along dimension p
			
		The sampling area is then enlarged a bit by \c enlargeAreaFactor to avoid that sampling area's extreme coincide with input MCMC samples.
		*/
		void createDefaultArea();
};
#endif