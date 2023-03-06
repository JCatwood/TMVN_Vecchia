#include <RcppEigen.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <cmath>
#include "mvphi.h"


using namespace std;
using namespace Eigen;


int mvndns(
    const Eigen::VectorXd &a, const Eigen::VectorXd &b, 
    const Eigen::MatrixXi &NN, const Eigen::MatrixXd &muCoeff, 
    const Eigen::VectorXd &condSd, const Eigen::VectorXd &beta, int NLevel1, int NLevel2)
{
	int n = a.rows();  // MVN dim
	int N = NLevel2 / 2;  
    int m = NN.cols() - 1;
	NLevel2 = N * 2;
	VectorXd p_L1(NLevel2);
	VectorXd psi_L2(NLevel1);
	MatrixXd MC_grid(n, N);
	VectorXd MC_rnd(n);
	MatrixXd MC_samp(n, NLevel2);
	VectorXd prime(n);
    VectorXd a_bat(NLevel2);
    VectorXd b_bat(NLevel2);
    VectorXd pnorm_at_a(NLevel2);
    VectorXd pnorm_at_b(NLevel2);
    VectorXd pnorm_diff(NLevel2);
    VectorXd cdf_MC_samp(NLevel2);
    MatrixXd X(n, NLevel2);
    VectorXd mu(NLevel2);
    VectorXd lnNpr_sum(NLevel2);
    VectorXd inner_prod(NLevel2);
    Matrix<int, n, m> cond_ind = NN.block(0, 1, n, m) - 1;

	// generate n prime numbers
        primes(5*(n+1)*log((double)(n+1)+1)/4, n, prime);

    // MC_grid = sqrt(prime) * one2N
    transform(prime, prime + n, MC_grid.data(), [](int x){
        return sqrt((double) x);});
    for(int i = 1; i < N; i++)
        transform(MC_grid.data(), MC_grid.data()+n, MC_grid.col(i-1).data(),
            MC_grid.col(i).data(), [](double x1,double x2){return x1 + x2;});

    // Level 1 iteration
    for(int k = 0; k < NLevel1; k++){
        lnNpr_sum.setZeros();
        inner_prod.setZeros();
        // generate MC_rnd from R RNG
        GetRNGstate();
        for_each(MC_rnd.data(), MC_rnd.data()+n, [](double &x){
            x = unif_rand();});
        PutRNGstate();
        // Fill in first half of MC_samp
        for(int j = 0; j < N; j++)
            transform(
                MC_grid.col(j).data(), MC_grid.col(j).data()+n, 
                MC_rnd.data(), MC_samp.col(j).data(), 
                [](double &xFix, double &xr){return xFix + xr;});
        for_each(MC_samp.data(), MC_samp.data() + n * N, 
            [](double &x){x = abs(2.0 * (x - int(x)) - 1.0);});
        // Fill in second half of MC_samp
        transform(
            MC_samp.data(), MC_samp.data() + n * N, 
            MC_samp.data() + n * N, [](double &x){return 1 - x;});

        // Level 2 iteration
        for(int i = 0; i < n; i++){
            // Copy a and b by N times
            fill(a_bat.data(), a_bat.data() + NLevel2, a(i));
            fill(b_bat.data(), b_bat.data() + NLevel2, b(i));
            // Compute mu
            if(i > 0){
                mu = muCoeff.row(i) * X(cond_ind.row(i), placeholders::all);
                a_bat.noalias() = a_bat - mu;
                b_bat.noalias() = b_bat - mu;
            }
            // scale a_bat and b_bat
            a_bat.noalias() = a_bat / condSd(i);
            b_bat.noalias() = b_bat / condSd(i);
            // substract beta
            a_bat.noalias() = a_bat - beta(i);
            b_bat.noalias() = b_bat - beta(i);
            // compute 1d normal cdf
            lc_vdCdfNorm(NLevel2, a_bat.data(), pnorm_at_a.data());
            lc_vdCdfNorm(NLevel2, b_bat.data(), pnorm_at_b.data());
            pnorm_diff = pnorm_at_b - pnorm_at_a;
            // sample X
            cdf_MC_samp = pnorm_at_a + MC_samp.row(i).cwiseProduct(pnorm_diff);
            lc_vdCdfNormInv(NLevel2, cdf_MC_samp.data(), 1, X.row(i).data(), n);
            X.row(i) = X.row(i) * condSd(i) + mu + beta(i) * condSd(i);
            // compute the i-th summand of psi
            lnNpr_sum.noalias() = lnNpr_sum + pnorm_diff.log();
            inner_prod.noalias() = inner_prod + (X.row(i) - mu) * beta(i) / condSd(i);
        }

    }
}
