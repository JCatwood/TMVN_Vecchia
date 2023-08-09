// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "mvphi.h"

#include <iostream>


using namespace Rcpp;
using namespace std;


void update_cond_pars(const arma::mat &corrMatSub, const arma::vec &corrVecSub,
	arma::mat &mean_coeff, arma::vec &cond_var, int j)
{
	int m = corrVecSub.n_elem;
	arma::mat mat_inv = arma::inv_sympd(corrMatSub);
	mean_coeff(j, arma::span(0, m - 1)) = (mat_inv * corrVecSub).t();
	cond_var(j) = 1.0 - arma::dot(mean_coeff(j, arma::span(0, m - 1)), 
		corrVecSub);
}

// [[Rcpp::export]]
List univar_order_vecc(arma::vec a, arma::vec b, arma::mat corrMat, int m)
{
	int n = a.n_elem;
	arma::umat NN(n, m, arma::fill::none);
	arma::mat mean_coeff(n, m, arma::fill::zeros);
	arma::mat dist(n, m, arma::fill::ones);
	arma::vec pnorm_at_a(n);
	arma::vec pnorm_at_b(n);
	arma::vec pnorm_diff(n);
	arma::vec cond_expectation(n);
	arma::vec cond_mean(n, arma::fill::zeros);
	arma::vec cond_var(corrMat.diag());
	arma::uvec odr(n);
	for(int i = 0; i < n; i++) 
		odr(i) = i;
	for(int i = 0; i < n; i++){
		if(i > 0){
			for(unsigned int j = i; j < n; j++){
				if(i > m){
					if(dist(j, m - 1) > - corrMat(odr(j), odr(i - 1))){
						// update dist, NN, mean_coeff, and cond_var
						dist(j, m - 1) = - corrMat(odr(j), odr(i - 1));
						NN(j, m - 1) = i - 1;
						arma::uvec odr_row_j = arma::sort_index(dist(j, arma::span(0, m - 1)));
						NN(j, arma::span(0, m - 1)) = NN(arma::uvec({j}), odr_row_j);
						dist(j, arma::span(0, m - 1)) = dist(arma::uvec({j}),
							odr_row_j);
						arma::mat corr_mat_sub = corrMat(odr(NN(j, arma::span(0, m - 1))), 
							odr(NN(j, arma::span(0, m - 1))));
						arma::vec corr_vec_sub = corrMat(arma::uvec({odr(j)}),
							odr(NN(j, arma::span(0, m - 1)))).t();
						update_cond_pars(corr_mat_sub, corr_vec_sub, mean_coeff, 
							cond_var, j);
						// compute cond_mean
						cond_mean(j) = 0;
						for(int k = 0; k < m; k++)
							cond_mean(j) += mean_coeff(j, k) * cond_expectation(NN(j, k));
					}
				}else{
					// update dist, NN, mean_coeff, and cond_var
					dist(j, i - 1) = - corrMat(odr(j), odr(i - 1));
					NN(j, i - 1) = i - 1;
					arma::uvec odr_row_j = arma::sort_index(dist(j, arma::span(0, i - 1)));
					NN(j, arma::span(0, i - 1)) = NN(arma::uvec({j}), odr_row_j);
					dist(j, arma::span(0, i - 1)) = dist(arma::uvec({j}), 
						odr_row_j);
					arma::mat corr_mat_sub = corrMat(odr(NN(j, arma::span(0, i - 1))), 
						odr(NN(j, arma::span(0, i - 1))));
					arma::vec corr_vec_sub = corrMat(arma::uvec({odr(j)}),
						odr(NN(j, arma::span(0, i - 1)))).t();
					update_cond_pars(corr_mat_sub, corr_vec_sub, mean_coeff, 
						cond_var, j);
					// compute cond_mean
					cond_mean(j) = 0;
					for(int k = 0; k < i; k++)
						cond_mean(j) += mean_coeff(j, k) * cond_expectation(NN(j, k));
				}
			}
		}
		arma::vec a_sub = a.subvec(i, n - 1) - cond_mean.subvec(i, n - 1);
		arma::vec b_sub = b.subvec(i, n - 1) - cond_mean.subvec(i, n - 1);
		a_sub /= sqrt(cond_var.subvec(i, n - 1));
		b_sub /= sqrt(cond_var.subvec(i, n - 1));
		lc_vdCdfNorm(n - i, a_sub.memptr(), pnorm_at_a.memptr() + i);
	    lc_vdCdfNorm(n - i, b_sub.memptr(), pnorm_at_b.memptr() + i);
	    arma::vec pnorm_diff = pnorm_at_b.subvec(i, n - 1) - 
	    	pnorm_at_a.subvec(i, n - 1);
	    int j_hat = pnorm_diff.index_min();
	    int i_hat = j_hat + i;
	    // compute cond_expectation
	    double dnorm_diff = (exp(- b_sub(j_hat) * b_sub(j_hat) / 2) - 
	    	exp(- a_sub(j_hat) * a_sub(j_hat) / 2)) / sqrt(2 * M_PI);
	    if(pnorm_diff(j_hat) < 1e-20) {
	    	if(isfinite(a_sub(j_hat))) {
	    		cond_expectation(i) = a_sub(j_hat) * 
	    			sqrt(cond_var(i_hat)) + cond_mean(i_hat);
	    	} else if(isfinite(b_sub(j_hat))) {
	    		cond_expectation(i) = b_sub(j_hat) * 
	    			sqrt(cond_var(i_hat)) + cond_mean(i_hat);
	    	} else {
	    		throw invalid_argument(
	    			"a[i_hat] and b[i_hat] are equal and both infinite\n");
	    	}
	    } else {
	    	cond_expectation(i) = - dnorm_diff / pnorm_diff(j_hat) * 
	    		sqrt(cond_var(i_hat)) + cond_mean(i_hat);
	    }
	    
	    // switch pairs in a, b, odr
	    swap(a(i), a(i_hat));
	    swap(b(i), b(i_hat));
	    swap(odr(i), odr(i_hat));
	    NN.swap_rows(i, i_hat);
	    mean_coeff.swap_rows(i, i_hat);
	    dist.swap_rows(i, i_hat);
	    swap(cond_mean(i), cond_mean(i_hat));
	    swap(cond_var(i), cond_var(i_hat));
	}
	return List::create(Named("order") = odr, Named("nn") = NN,
		Named("cond_mean_coeff") = mean_coeff, Named("cond_var") = cond_var);
}
