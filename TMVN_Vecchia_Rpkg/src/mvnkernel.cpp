#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include "mvphi.h"


using namespace Eigen;
using namespace std;


void primes(int n, int sz, int *primeVec)
{
    int idx = 0;
    if(n > 2 && sz > 0)
    {
        primeVec[idx] = 2;
        idx++;
        if(idx == sz)
            return;
        for(int i = 3; i <= n; i++)
        {
            int sqroot = sqrt(i);
            bool prime = true;
            for(int j = 0; j < idx; j++)
            {
                if(primeVec[j] > sqroot)
                    break;
                if(i % primeVec[j] == 0)
                {
                    prime = false;
                    break;
                }
            }
            if(prime)
            {
                primeVec[idx] = i;
                idx++;
                if(idx == sz)
                    return;
            }
        }   
    }
}


// [[Rcpp::export]]
Eigen::VectorXd mvndns(
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
	VectorXi prime(n);
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
    MatrixXi cond_ind = NN.block(0, 1, n, m).array() - 1;

	// generate n prime numbers
    primes(5*(n+1)*log((double)(n+1)+1)/4, n, prime.data());

    // MC_grid = sqrt(prime) * one2N
    transform(prime.data(), prime.data() + n, MC_grid.data(), [](int x){
        return sqrt((double) x);});
    for(int i = 1; i < N; i++)
        transform(MC_grid.data(), MC_grid.data()+n, MC_grid.col(i-1).data(),
            MC_grid.col(i).data(), [](double x1,double x2){return x1 + x2;});

    // Level 1 iteration
    for(int k = 0; k < NLevel1; k++){
        lnNpr_sum.setZero();
        inner_prod.setZero();
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
                mu.setZero();
                for(int j = 0; j < m; j++){
                    if(muCoeff(i, j) != 0){
                        mu.noalias() = mu + muCoeff(i, j) * X.row(cond_ind(i, j));
                    }
                }
                a_bat.noalias() = a_bat - mu;
                b_bat.noalias() = b_bat - mu;
            }
            // scale a_bat and b_bat
            a_bat.noalias() = a_bat / condSd(i);
            b_bat.noalias() = b_bat / condSd(i);
            // substract beta
            a_bat.array() -= beta(i);
            b_bat.array() -= beta(i);
            // compute 1d normal cdf
            lc_vdCdfNorm(NLevel2, a_bat.data(), pnorm_at_a.data());
            lc_vdCdfNorm(NLevel2, b_bat.data(), pnorm_at_b.data());
            pnorm_diff = pnorm_at_b - pnorm_at_a;
            // sample X
            cdf_MC_samp = pnorm_at_a + MC_samp.row(i).cwiseProduct(pnorm_diff);
            lc_vdCdfNormInv(NLevel2, cdf_MC_samp.data(), 1, X.row(i).data(), n);
            X.row(i).noalias() = X.row(i) * condSd(i);
            X.row(i).noalias() = X.row(i) + mu;
            X.row(i).array() += beta(i) * condSd(i);
            // compute the i-th summand of psi
            lnNpr_sum.noalias() = lnNpr_sum + pnorm_diff.log();
            inner_prod.noalias() = inner_prod + (X.row(i) - mu) * beta(i) / condSd(i);
        }
        psi_L2 = - inner_prod.array() + lnNpr_sum.array() + 0.5 * beta.squaredNorm();
        p_L1(k) = psi_L2.exp().mean();
    }
    return p_L1;
}
