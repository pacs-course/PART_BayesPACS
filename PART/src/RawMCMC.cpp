#include <stdexcept>
#include <cassert>
#include <cmath>
#include "RawMCMC.hpp"

using namespace arma;

// Constructor from vector of pointers to matrices of MCMC samples drawn on each subset
RawMCMC::RawMCMC(const std::vector<const mat*> & subchains):m_M(subchains.size()){
	
	if(m_M<1){
		throw std::length_error ("Too few subsets!!");	
	}
	
	const uword d=(subchains[0])->n_cols; // number of features
	if(d<1){
		throw std::length_error ("Too few features!!");	
	}
	
	uword nMCsamples=0; // number of samples to be stored in m_samples
	for(const mat * const i: subchains){
		
		if(i->n_rows<1){
			throw std::length_error ("Too few rows in MCMC matrix!!");	
		}
		
		nMCsamples += i->n_rows;
		
	}
	
	// set m_sample, m_mark and m_sumLogN
	m_samples.set_size(nMCsamples,d);
	m_mark.set_size(nMCsamples);
	m_sumLogN = 0.0;
	uword m = 0;
	uword id_f = 0;
	for(const mat * const i: subchains){
		const uword n = i->n_rows; // number of rows of current MCMC samples' matrix
		m_sumLogN += std::log(n);
		const uvec tmp_rows = linspace<uvec>(id_f,id_f+n-1,n); // indices of m_samples' and m_mark's rows that need to be filled
		m_samples.rows(tmp_rows)=(*i); // fill m_samples with MCMC samples of interest
		m_mark(tmp_rows)=m*ones<uvec>(n); // fill m_mark with index of subset from which MCMC samples have been drawn
		m++; // move to next subset index
		id_f = id_f + n;
	}
	
	// create sampling area m_area
	createDefaultArea();
	
};

// get method for m_M attribute
const uword RawMCMC::get_n_subsets() const { return m_M; };

// get method for m_sumLogN attribute
const double RawMCMC::get_sumLogN() const { return m_sumLogN; };

// get method for m_samples attribute
const mat & RawMCMC::get_samples() const { return m_samples; };

// get method for m_mark attribute
const uvec & RawMCMC::get_mark() const { return m_mark; };

// get method for m_area attribute
const mat & RawMCMC::get_area() const { return m_area; };

// set method for m_M attribute
void RawMCMC::set_n_subsets(const uword nSubsets) { m_M = nSubsets; };

// set method for m_sumLogN attribute
void RawMCMC::set_sumLogN(const double sum_log_N) { m_sumLogN = sum_log_N; };
	
// set-method for m_samples attribute
void RawMCMC::set_samples(const mat & samples) { m_samples = samples; };
	
// set-method for m_mark attribute
void RawMCMC::set_mark(const uvec & mark) { m_mark = mark; };
	
// set-method for m_area attribute
void RawMCMC::set_area(const mat & area) { m_area = area; };

// Method to compute and set sampling area from MCMC samples
void RawMCMC::createDefaultArea(){
	
	assert(enlargeAreaFactor>0.0 && "Error in default area: enlargeAreaFactor must be a positive number!!");
	
	rowvec row_min = min(m_samples,0), row_max = max(m_samples,0); // min and max of MCMC samples in m_samples for each column
	assert(all(row_max - row_min>=0.0) && "Error in default area: left extreme cannot be higher than right extreme!!");
	
	// enlarge sampling area by enlargeAreaFactor
	const uword d = m_samples.n_cols;
	for(uword p=0; p<d; p++){

		if(row_min(p)>0.0){
			row_min(p) /= enlargeAreaFactor;
		}else{
			row_min(p) *= enlargeAreaFactor;
		}

		if(row_max(p)>0.0){
			row_max(p) *= enlargeAreaFactor;
		}else{
			row_max(p) /= enlargeAreaFactor;
		}
		
	}
	assert(all(row_max - row_min>=0.0) && "Error in default area: left extreme cannot be higher than right extreme!!");
	
	// store left and rigth boundaries into sampling area
	m_area.set_size(2,d);
	m_area.row(0)=row_min;
	m_area.row(1)=row_max;
	assert( all(row_max - max(m_samples,0)>=0.0) && all(min(m_samples,0)-row_min>=0.0) && "Error in default area: left or right extreme of default area do not include all MCMC samples!!");

};