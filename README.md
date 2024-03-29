README
================

# its2es

This package implements interrupted time series analysis for both
continuous and count outcomes, and quantifies the associated effect
size, as described in Effect size quantification for interrupted time
series analysis: Implementation in R and analysis for Covid-19 research.
The main functions fit an ITS regression model, and then use the fitted
values and the model-based counterfactual values to quantify the effect
size (Cohen’s d for continuous outcomes and relative risk for count
outcomes). An example describing how to install and use this package is
described below. A more detailed tutorial, including the data analysis
described in the paper, is also available with this package (Rmd + pdf
file).

## Installation

You can install the package from its [GitHub
repository](https://github.com/Yael-Travis-Lumer/its2es/). You first
need to install the [devtools](https://github.com/r-lib/devtools)
package.

``` r
install.packages("devtools",repos = "http://cran.us.r-project.org")
```

Then install its2es using the `install_github` function in the
[devtools](https://github.com/r-lib/devtools) package.

``` r
library(devtools)
install_github("Yael-Travis-Lumer/its2es")
```

## Example

1.  Load library and Israel all-cause mortality data (discussed in
    paper)

``` r
library(its2es)
data <- Israel_mortality
```

2.  Define formula and intervention start index for the Covid-19 period

``` r
form <- as.formula("percent ~ time")
intervention_start_ind <- which(data$Year==2020 & data$Month==3)
```

3.  Fit a linear regression ITS model to the mortality percent

``` r
fit <- its_lm(data=data,form=form,time_name = "time",intervention_start_ind=intervention_start_ind, freq=12,seasonality= "full", impact_model = "full",counterfactual = TRUE, print_summary = FALSE)
```

    ## Cohen's d   2.5% CI  97.5% CI   P-value 
    ##  1.038391  0.332192  1.715101  0.002500

4.  Plot predicted values and counterfactual values

``` r
p <- plot_its_lm(data=fit$data,intervention_start_ind=intervention_start_ind, y_lab="All-cause mortality percent", response="percent", date_name= "Date")
p
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
