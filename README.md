# MLDID

The **MLDID** package computes average and conditional average treatment effects on the treated (ATT and CATT), allowing for them to be used to study detailed drivers of treatment effect heterogeneity, see:

Hatamyar, J., Huber, M., Rocha, R. Kreif, N. *Machine Learning for Staggered Difference-in-Differences and Dynamic Treatment Effect Heterogeneity* (2023) [Arxiv](https://arxiv.org/abs/2310.11962)

## Installation 

You can install the package using 

```
# install.packages("devtools")
devtools::install_github("jhatamyar/MLDID")
```

## Basics 

This code demonstrates the basic functionality of the package using the data on minimum wage and county-level employment from Callaway & Sant'Anna (2021). Note that there is no treatment effect heterogeneity in this data, and it is used here for demonstrative purposes. Vignette for more detailed usage forthcoming. 

``` R
library(MLDID)
library(did) ## to get the data
data(mpdta)
```
