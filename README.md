# its2es
This package implements interrupted time series analysis for both continuous and count outcomes, and quantifies the associated effect size, as described in Effect Size Quantification for Interrupted Time Series Analysis: Implementation in R for Covid-19 Research.     The main functions fit an ITS regression model, and then use the fitted values and the model-based counterfactual values to quantify the effect size (Cohen's d for continuous outcomes and relative risk for count outcomes).

## Installation
You can install the package from its [GitHub repository](https://github.com/Yael-Travis-Lumer/its2es/). You first need to install the [devtools](https://github.com/r-lib/devtools) package.
```{r}
install.packages("devtools")
```
Then install its2es using the `install_github` function in the [devtools](https://github.com/r-lib/devtools) package.
```{r}
library(devtools)
install_github("Yael-Travis-Lumer/its2es")
```

## Example
1. Load library and unemployment data (discussed in paper)
```{r}
library(its2es)
data <- unemployed
```

2. Define formula and intervention start index for the Covid-19 period
```{r}
form <- as.formula("percent ~ time")
intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
```
3. Fit a linear regression ITS model
```{r}
fit <- its_lm(data=data,form=form,time_name = "time",intervention_start_ind=intervention_start_ind, freq=12,seasonality= "none", impact_model = "full",counterfactual = TRUE)
```
4. Plot predicted values and counterfactual values
```{r}
p <- plot_its_lm(data=fit$data,intervention_start_ind,y_lab="Unemployment percent",response="percent", date_name= "dt")
```

