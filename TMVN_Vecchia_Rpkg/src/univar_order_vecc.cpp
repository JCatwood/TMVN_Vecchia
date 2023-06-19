#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "mvphi.h"


using namespace Rcpp;
using namespace std;


void update_cond_pars(const arma::mat &corrMatSub, const arma::uvec &corrVecSub,
	arma::mat &mean_coeff, arma::vec &cond_var, int j)
{

}


List univar_order_vecc(arma::vec a, arma::vec b, arma::mat corrMat, int m)
{
	int n = a.length();
	arma::umat NN(n, m, arma::fill::none);
	arma::mat mean_coeff(n, m, arma::fill::zeros);
	arma::mat dist(n, m, arma::fill::zeros);
	arma::vec pnorm_at_a(n);
	arma::vec pnorm_at_b(n);
	arma::vec pnorm_diff(n);
	arma::vec cond_expectation(n);
	arma::vec cond_mean(n, arma::fill::zeros);
	arma::vec cond_var(covMat.diag());
	arma::uvec odr(n);
	for(int i = 0; i < n; i++) 
		odr(i) = i;
	for(int i = 0; i < n; i++){
		if(i > 0){
			if(i > m){
				// update NN, mean_coeff, and cond_var
				
				// compute cond_mean
			}else{
				for(int j = i; j < n; j++){
					dist(j, i - 1) = - corrMat(odr(j), odr(i - 1));
					NN.subvec(0, i - 1) = sort_index(dist.subvec(0, i - 1));
					arma::mat corr_mat_sub = corrMat(odr(NN.subvec(0, i - 1)), 
						odr(NN.subvec(0, i - 1)));
					arma::vec corr_vec_sub = corrMat(odr(j), 
						odr(NN.subvec(0, i - 1)));
					update_cond_pars(corr_mat_sub, corr_vec_sub, mean_coeff, 
						cond_var, j);
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

        // switch pairs in a, b, odr
        

	}
}
