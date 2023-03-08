#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <numeric>
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



/*
    For 2d arrays, row number correspond to MVN dim and col number 
        correspond to MC sample. 
    Indexing is row-major for 2d arrays.
*/
// [[Rcpp::export]]
NumericVector mvndns(
    const NumericVector &a, const NumericVector &b, 
    const IntegerMatrix &NN, const NumericMatrix &muCoeff, 
    const NumericVector &condSd, const NumericVector &beta, 
    int NLevel1, int NLevel2)
{
    int n = a.rows();  // MVN dim
    int N = NLevel2 / 2;  
    int m = NN.cols() - 1;
    NLevel2 = N * 2;
    NumericVector p_L1(NLevel1);
    double * psi_L2 = new double[NLevel2];
    double * MC_grid = new double[n * N];
    double * MC_rnd = new double[n];
    double * MC_samp = new double[n * NLevel2];
    int * prime = new int[n];
    double * a_bat = new double[NLevel2];
    double * b_bat = new double[NLevel2];
    double * pnorm_at_a = new double[NLevel2];
    double * pnorm_at_b = new double[NLevel2];
    double * pnorm_diff = new double[NLevel2];
    double * cdf_MC_samp = new double[NLevel2];
    double * X = new double[n * NLevel2];
    double * mu = new double[NLevel2];
    double * lnNpr_sum = new double[NLevel2];
    double * inner_prod = new double[NLevel2];
    int * cond_ind = new int[n * m]; 

    // copy nearest neighbors into cond_ind
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            cond_ind[i * m + j] = NN(i, j + 1) - 1;

    // generate n prime numbers
    primes(5*(n+1)*log((double)(n+1)+1)/4, n, prime);

    // MC_grid = sqrt(prime) * one2N
    for(int i = 0; i < n; i++)
        MC_grid[i * N] = sqrt((double) prime[i]);
    for(int i = 0; i < n; i++)
        for(int j = 1; j < N; j++)
            MC_grid[i * N + j] = MC_grid[i * N + j - 1] + MC_grid[i * N];

    // Level 1 iteration
    for(int k = 0; k < NLevel1; k++){
        fill(lnNpr_sum, lnNpr_sum + NLevel2, 0.0);
        fill(inner_prod, inner_prod + NLevel2, 0.0);
        // generate MC_rnd from R RNG
        GetRNGstate();
        for_each(MC_rnd, MC_rnd + n, [](double &x){x = unif_rand();});
        PutRNGstate();
        // Fill in MC_samp
        for(int i = 0; i < n; i++)
            for(int j = 0; j < N; j++){
                double tmp_val = MC_grid[i * N + j] + MC_rnd[i];
                MC_samp[i * NLevel2 + j] = 
                    abs(2.0 * (tmp_val - int(tmp_val)) - 1.0);
                MC_samp[i * NLevel2 + N + j] = 1.0 - MC_samp[i * NLevel2 + j];
            }

        // Level 2 iteration
        for(int i = 0; i < n; i++){
            // Copy a and b by NLevel2 times
            fill(a_bat, a_bat + NLevel2, a[i]);
            fill(b_bat, b_bat + NLevel2, b[i]);
            // Compute mu
            if(i > 0){
                fill(mu, mu + NLevel2, 0.0);
                for(int j = 0; j < m; j++){
                    int cond_ind_j = cond_ind[i * m + j];
                    double mu_coeff_j = muCoeff(i, j);
                    double * X_row_i = X + i * NLevel2;
                    if(mu_coeff_j != 0)
                        for(int j2 = 0; j2 < NLevel2; j2++)
                            mu[j2] += mu_coeff_j * X_row_i[j2];
                }
                for(int j = 0; j < NLevel2; j++){
                    a_bat[j] -= mu[j];
                    b_bat[j] -= mu[j];
                }
            }
            // update a_bat and b_bat
            double cond_sd_i = condSd[i];
            double beta_i = beta[i];
            for(int j = 0; j < NLevel2; j++){
                a_bat[j] /= cond_sd_i;
                b_bat[j] /= cond_sd_i;
                a_bat[j] -= beta_i;
                b_bat[j] -= beta_i;
            }
            // compute 1d normal cdf
            lc_vdCdfNorm(NLevel2, a_bat, pnorm_at_a);
            lc_vdCdfNorm(NLevel2, b_bat, pnorm_at_b);
            transform(pnorm_at_a, pnorm_at_a + NLevel2, pnorm_at_b, pnorm_diff,
                [](double &x, double &y){return y - x;});
            // sample X and compute i-th summand of psi
            double * MC_samp_row_i = MC_samp + i * NLevel2;
            for(int j = 0; j < NLevel2; j++){
                cdf_MC_samp[j] = pnorm_at_a[j] + MC_samp_row_i[j] * 
                    pnorm_diff[j];
            }
            double * X_row_i =  X + i * NLevel2;
            lc_vdCdfNormInv(NLevel2, cdf_MC_samp, X_row_i);
            for(int j = 0; j < NLevel2; j++){
                X_row_i[j] = X_row_i[j] * cond_sd_i + mu[j] + 
                    beta_i * cond_sd_i;
                lnNpr_sum[j] += log(pnorm_diff[j]);
                inner_prod[j] += (X_row_i[j] - mu[j]) * beta_i / cond_sd_i;
            }
        }
        double beta_sq_norm = inner_product(beta.begin(), beta.end(), 
            beta.begin(), 0.0);
        for(int j = 0; j < NLevel2; j++){
            psi_L2[j] = exp(- inner_prod[j] + lnNpr_sum[j] + 
                0.5 * beta_sq_norm);
        }
        p_L1[k] = accumulate(psi_L2, psi_L2 + NLevel2, 0.0) / NLevel2;
    }

    delete[] psi_L2;
    delete[] MC_grid;
    delete[] MC_rnd;
    delete[] MC_samp;
    delete[] prime;
    delete[] a_bat;
    delete[] b_bat;
    delete[] pnorm_at_a;
    delete[] pnorm_at_b;
    delete[] pnorm_diff;
    delete[] cdf_MC_samp;
    delete[] X;
    delete[] mu;
    delete[] lnNpr_sum;
    delete[] inner_prod;
    delete[] cond_ind;

    return p_L1;
}
