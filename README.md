# This is the Git directory for the paper **Linear-Cost Vecchia Approximation of Multivariate Normal Probabilities**

## The paper 
The corresponding paper can be found [here](https://arxiv.org/abs/2311.09426).

## VeccTMVN R package
The methods developed in the paper are implemented in the R package `VeccTMVN`. To install the package, simply run
```
install.packages("VeccTMVN")
```

## Contents
  - experiments: The codes for reproducing the results in the paper. All scripts can be run independently, i.e., there is no requirement on the order for running the R scripts in this directory. However, to run the scripts, one first needs to install the dependent R packages; please refer to the "dependency installation" Section.
    - `ordering_bias.R` for Figure 1 and Figure C.2
    - `Vecc_bias_lowdim.R` for Figure 2
    - `Vecc_bias_highdim.R` for Figure 3
    - `Vecc_bias_const_corr.R` for Figure 4
    - `TMVN_sim_low.R` for Figure 5
    - `TMVN_parm_est_high.R` for Figure 6
    - `PTMVN_sim_low.R` for Figure 7
    - `PTMVN_sim_high.R` for the RMSE reported in the 6400 dimensional example in Section 7.2 and D.3
    - `PCE_modeling.R` for Figure 8
    - `Vecc_bias_lowdim_T.R` for Figure B.1
    - `results` this folder stores results from running some of the above R scripts
    
    **Some scripts can take long. Alternatively, one can comment out the computation-heavy section as instructed by the comments inside each script. Then the scripts will use the previously computed results to obtain the figures.**

  - external_Rpkg: R packages needed for `experiments` that are not on CRAN. Specifically,
    - CDFApprox: code from the paper ''A Vecchia approximation for high-dimensional Gaussian cumulative distribution functions arising from spatial data'' (https://doi.org/10.1080/00949655.2021.2016759)
    - TruncatedNormalBeta: Modified version of the R package `TruncatedNormal` that implements Botev (2017). The modification removes the univariate reordering to demonstrate the impact of variable reordering as shown in Figure 1
  - funcs: some peripheral functions needed for running the experiments

## Dependency installation
Tested out on a Linux system (UH Sabine Cluster). Should work on other platforms but modification of installation command may be needed.

R code
```
install.packages(c("VeccTMVN","sf","spData","GpGp","scoringRules","fields",
  "RColorBrewer","mvtnorm","TruncatedNormal",
  "CensSpatial","ggplot2","tidyr","scales","tlrmvnmvt","dplyr"))
```

Bash script
```
cd external_Rpkg/
R CMD INSTALL CDFApprox/
R CMD INSTALL TruncatedNormalBeta/
```

## Data
The results mostly use simulated data. There is only exception that the script `PCE_modeling.R` uses a real dataset, the Tetrachloroethylene concentration dataset from the United States Geological Survey (USGS), which is publically available and provided in this repository at `experiments/PCE.RData` for convenience.
	
## A simple example of VeccTMVN
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
