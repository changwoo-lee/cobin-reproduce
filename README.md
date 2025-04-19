
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cobin-reproduce

<!-- badges: start -->
<!-- badges: end -->

Code to reproduce figures, simulation results, and case study results
from the paper

> Anonymous authors (2025). Scalable and robust regression models for
> continuous proportional data. \[submitted\]

## Prerequisites

\[to be added very soon after anonymizing\] The R package `cobin`
contains all necessary functions for cobin and micobin regression
models. Download `cobin_1.0.0.0.tar.gz` and install cobin R package from
the source using the following command in R:

``` r
# install.packages("devtools")
install.packages("path_to_downloaded_file/cobin_1.0.0.0.tar.gz", repos = NULL, type="source")
```

Test if the package is installed correctly, e.g. Kolmogorov-Gamma
sampler:

``` r
library(microbenchmark)
library(BayesLogit)
library(cobin)
nsample = 1e6 # million sample
microbenchmark(
  BayesLogit::rpg(nsample, 1, 2), # Polya-Gamma sampler
  cobin::rkgcpp(nsample, 1, 2) # Kolmogorov-Gamma sampler
)
#> Unit: milliseconds
#>                            expr      min       lq     mean   median       uq
#>  BayesLogit::rpg(nsample, 1, 2) 208.4826 257.4689 263.7239 260.7093 267.7383
#>    cobin::rkgcpp(nsample, 1, 2) 509.6457 517.1453 525.8275 520.9337 529.2793
#>       max neval
#>  337.4130   100
#>  594.6312   100
```

## Codes for reproducing simulation results

### Simulation 1: consistency under model misspecification

- `Sec5_simul1/betareg_cobit.R`: code for finding beta regression MLE
  with cobit link; derived from betareg package
- `Sec5_simul1/datagen_consistency.R`: code for synthetic data
  generation.
- `Sec5_simul1/runsim_consistency.R`: code for running simulation.

### Simulation 2: scalablity and robustness using spatial regression under model misspecification

There are total 4 scenarios: rho01_n250, rho01_n500, rho02_n250,
rho02_n500

- `Sec5_simul2/betareg_spatial.stan`: code for spatial beta regression
  model with cobit link in Stan
- `Sec5_simul2/datagen_spatial.R`: code for synthetic data generation.
- `Sec5_simul2/[scenario]/runbeta_[scenario]`: code for running spatial
  beta regression.
- `Sec5_simul2/[scenario]/runcobin_[scenario]`: code for running spatial
  cobin regression.
- `Sec5_simul2/[scenario]/runmicobin_[scenario]`: code for running
  spatial micobin regression.

## Codes for reproducing benthic macroinvertibrate multimetric index (MMI) analysis

### Data preparation, exploratory data analysis

- `Sec6_mmicasestudy/mmi_lakecat_csv`: MMI and lake watershed covariate
  for 950 lakes
- `Sec6_mmicasestudy/lakecat-over40000m2.csv`: lake watershed covariate
  for 50k+ lakes for prediction
  - See dataprocess_readme.txt for data collecting procedure from EPA
    website
- `Sec6_mmicasestudy/fig_mmidata.R`: code for reproducing figure 1 (MMI
  data and urban cover)

### Data analysis and prediction

- `Sec6_mmicasestudy/results_main_n949`: folder containing MMI data
  analysis of 949 lakes, excluding one lake with MMI = 0
- `Sec6_mmicasestudy/results_supp_n947`: folder containing MMI data
  analysis of 947 lakes, excluding two additional lakes from 949 lakes
  for sensitivity analysis
- `Sec6_mmicasestudy/results_supp_n950`: folder containing MMI data
  analysis of 950 lakes (only micobin regression) for sensitivity
  analysis
- `Sec6_mmicasestudy/fig_prediction.R`: code for reproducing figure 3
  (prediction)
- `Sec6_mmicasestudy/fig_qresiduals.R`: code for reproducing figure 4
  (quantile residual)

## Codes for reproducing other figures

- Figure 2 (density of cobin and micobin):
  - `others/density.R`
- Figure 2 (link function and variance function comparison):
  - `others/linkft.R`
