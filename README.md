## CLRtools: Diagnostic Tools for Logistic and Conditional Logistic Regression

This repository contains the source code for the CLRtools R package. CLRtools provides functions to support the structured development of logistic regression and conditional logistic regression models. The package also includes diagnostic tools for both modelling frameworks, including residual based diagnostic plots, influence diagnostics, and Bayesian inspired diagnostic plots to support model assessment.

The package includes vignettes that demonstrate the complete modelling workflow. The logistic regression vignette uses the GLOW500 dataset, which examines risk factors for fracture in women over 50. The conditional logistic regression vignette uses the GLOW11M dataset, a 1:1 matched case control dataset derived from GLOW500, where each woman with a fracture is matched to a woman of the same age without a fracture.

Contributions and suggestions are welcome.

## Installation and loading

- Install from [CRAN](https://cran.rstudio.com/web/packages/CLRtools/) as
  follow:

``` r
install.packages("ggpubr")
```

- Or, install the latest version from
  [GitHub](https://github.com/brendacontla/CLRtools) as follow:

``` r
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("brendacontla/CLRtools")
```
