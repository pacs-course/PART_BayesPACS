#include <stdexcept>
#include "computeMeanCov.hpp"
using namespace arma;

// Compute mean and covariance of MCMC samples
void computeMeanCov (const mat & samples, const uvec & mark, const uword M, rowvec & mu, mat & Sigma){
	
	if(samples.n_rows<1 || mark.n_rows<1){
		throw std::length_error ("Too few rows in MCMC matrix!!");
	}
	
	if(samples.n_rows!=mark.n_rows){
		throw std::length_error ("Samples' and mark's dimension mismatch!!");
	}
	
	if(M<1){
		throw std::length_error ("Too few subsets!!!");
	}
	
	const uword d = samples.n_cols;
	if(d<1){
		throw std::length_error ("Too few features!!!");
	}

	// set sizes of covariance and mean
	Sigma.set_size(d,d);
	mu.set_size(d);
	
	if(mark.n_rows==1){ // mean and covariance undefined
		
		Sigma.fill(datum::nan);
		mu.fill(datum::nan);
	
	}else{
		
		uvec n_block(M); // n_block(m) = number of samples belonging to m-th subset
		for(uword m=0; m<M; m++){
			const uvec tmp = find(mark==m);
			n_block(m) = tmp.n_rows;
		}
		
		if(M==1 || !(all(n_block>d))){ // if samples drawn from only one subset OR at least one of the current subsets has less than d samples in it, use only simple estimates for mean and covariance of samples
			
			Sigma = cov(samples)/M;
			mu = mean(samples,0);
			
		}else{ // if there is more than one subset AND, for all subsets, the number of samples is higher than d, try block estimates
			
			mat aggInvSigma; // inverse of candidate covariance matrix
			aggInvSigma.zeros(d,d); 
			rowvec aggMu; // candidate mean vector
			aggMu.zeros(d);
			bool is_invertible = true;
			for(uword m=0; m<M; m++){
				
				const mat & sub_samples = samples.rows(find(mark==m)); // select samples that belong to subset m
				
				mat sub_invCov;
				is_invertible = inv(sub_invCov,cov(sub_samples)); // compute inverse of covariance of sub_samples
				if(!is_invertible){ // if covariance is not invertible avoid further checking
					break; 
				}
				
				aggInvSigma += sub_invCov;
				aggMu += mean(sub_samples,0)*sub_invCov;
				
			}
			
			mat cov_mat;
			vec eigval;
			mat eigvec;			
			if(is_invertible && inv(cov_mat,aggInvSigma) && eig_sym(eigval, eigvec, cov_mat) && all(eigval>=0.0) ){ // all sub covariances are invertible AND candidate covariance matrix is invertible and positive semi-definite
				
				Sigma = cov_mat;
				mu = aggMu*Sigma;
				
			}else{ // go back to simple estimates
				
				Sigma = cov(samples)/M;
				mu = mean(samples,0);
				
			}
			
		}
		
	}
	
};