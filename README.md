# This is the Git directory for the paper **Linear-Cost Vecchia Approximation of Multivariate Normal Probabilities**

## The paper 
The corresponding paper can be found [here](https://arxiv.org/abs/2311.09426).

## Installation
Simply run
```
install.packages("VeccTMVN")
```

## Directory Guidance

  - experiments: The codes for reproducing the results in the paper. The codes should be run after the VeccTMVN package is installed. Certain R packages on CRAN are needed. Packages in `external_Rpkg` are also needed.
  - external_Rpkg: R packages needed for `experiments` that are not on CRAN. Specifically,
    - CDFApprox: code from the paper ''A Vecchia approximation for high-dimensional Gaussian cumulative distribution functions arising from spatial data'' (https://doi.org/10.1080/00949655.2021.2016759)
    - TruncatedNormalBeta: Modified version of the R package `TruncatedNormal` that implements Botev (2017). The modification removes the univariate reordering 
  - funcs: some peripheral functions needed for running the experiments
	
## A simple example
```
library(GpGp)
library(VeccTMVN)
set.seed(123)
n1 <- 10
n2 <- 10
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- GpGp::matern15_isotropic(covparms, locs)
a <- rep(-Inf, n)
b <- -runif(n) * 2
m <- 30
est_Vecc <- VeccTMVN::pmvn(
  a, b, 0, locs,
  covName = "matern15_isotropic",
  covParms = covparms, m = m, verbose = F
)
```
