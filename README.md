# MLDID

The **MLDID** package computes average and conditional average treatment effects on the treated (ATT and CATT), allowing for them to be used to study detailed drivers of treatment effect heterogeneity, see:

Hatamyar, J., Kreif, N., Rocha, R., Huber, M. *Machine Learning for Staggered Difference-in-Differences and Dynamic Treatment Effect Heterogeneity* (2023) [Arxiv](https://arxiv.org/abs/2310.11962)

## Installation 

You can install the package using 

```
# install.packages("devtools")
devtools::install_github("jhatamyar/MLDID")
```

## Basics 

This code demonstrates the basic functionality of the package using the data on minimum wage and county-level employment from Callaway & Sant'Anna (2021). Note that it is used here for demonstrative purposes only. Vignette for more detailed usage forthcoming. 

``` R
library(MLDID)
library(did) ## to get the data
data(mpdta)
data <- mpdta
### There is only one covariate in the data so generate a noisy extra, as MLDID requires more than one covariate
data$noise <- rnorm(nrow(mpdta))
```

### The Group-Time ATT and CATT estimator 

The function `MLDID()` performs the main work of the package, implementing Algorithm 1 in the paper until the aggregation step. Functionality to estimate the nuisance functions using either rlearner or causal forest has been added, with the option to use the Superlearner to estimate delta, and the `tune_penalty` parameter can be set to `TRUE` to implement cross-validated parameter tuning for causal forest and the rlasso. 

```R
att_gt.mp <- MLDID(outcome = 'lemp',
                group = 'first.treat',
                time = 'year',
                id_name = 'countyreal',
                data,
                xformla = ~lpop + noise,
                tune_penalty = F,
                #nu_model = "cf",
                #sigma_model = "cf",
                #delta_model = "SL",
                )
```
The `att_gt.mp` object holds the estimated ATT and CATT for each group-time, as well as other information needed in post-processing steps. 

### Aggregating the estimates to Dynamic (time-to-event) time

The functions `dynamic_att()` and `dynamic_cates()` perform the aggregation steps in Equation 20 and 21 of the paper. The `dynamic_cates()` function can also aggregate the estimated robust scores.

```R
## aggregate the ATTs
ATT.dynamic.MP <- dynamic_attgt(att_gt.mp)
## print the dynamic estimates and the DRDID as in Callaway & Sant'Anna (2021) version
ATT.dynamic.MP[["dynamic.att.e"]]
ATT.dynamic.MP[["dynamic.att.e.csa"]]


## aggregate the CATTs and scores
cates.dynamic.MP <- dynamic_cates(att_gt.mp, type = "cates")
scores.dynamic.MP <- dynamic_cates(att_gt.mp, type = "scores")
```

### Heterogeneity Analysis 

The functions `BLP_eventtimes(), CLAN_glhtest() and CLAN_ttest()` allow for inference on the estimated CATTs and scores. The dynamic CATTs/scores must first be appended to the original data according to event-time using the function `het_prep()`. The function `BLP_summary()` provides a readable output of the coefficients for each time period and their standard errors. 

```R
## prepare the data for analysis. Note depending on the size of your data this may also be slow 
het.data.MP <- het_prep(att_gt.mp, cates.dynamic.MP)

## create a list of variables to test for heterogeneity, must be string
affected.mp <- c("lpop", "noise")
## formulate them for the input - note that here you could also just use "lpop + noise" as the argument for the functions instead of "affected_str.mp"
affected_str.mp <- paste(affected.mp, collapse = " + ")

## Run the BLP regressions: if error, reduce nperiods. 
BLP.bye.MP <- BLP_eventtimes(data = het.data.MP, nperiods = 3, rhs_formula = affected_str.mp)

## summarize the BLP output
BLPsummary(BLP.bye.MP, affected.mp)

## you can also check the summary for regression output at each event-time e using BLP.bye.MP[[e+1]][["coefficients"]]:
BLP.bye.MP[[1]][["coefficients"]] ## will be output for e=0
BLP.bye.MP[[2]][["coefficients"]] ## output for e=1 

## Run the CLANs two ways, first using the glh test as in Chernozhukov et al (2018), then a simple ttest of means of the most/least affected groups
CLAN.glh.cates.MP <- CLAN_glhtest(het.data.MP, affected=affected.mp)
CLAN.ttest.cates.MP <- CLAN_ttest(het.data.MP, affected=affected.mp)

## note there is no post-processing function for CLANS as of this update, but will be forthcoming
## can check the output by printing the object. Each matrix represents an event time. 
CLAN.glh.cates.MP 

## You can also visualize the BLP coefficients - we recommend interpreting with caution, as ideally the lpop variable should be discretized:
plotBLP(BLP.bye.MP, affected.mp)
```





![packagedemo](https://github.com/jhatamyar/MLDID/assets/31328293/8e012ac9-7dbf-4da9-9d2d-369988d93423)


### References
Callaway, B. and P. H. Sant’Anna (2021). *Difference-in-differences with multiple time periods*. In: Journal of Econometrics 225.2, pp. 200–230.

Chernozhukov, V., M. Demirer, E. Duflo, and I. Fernandez-Val (2018). *Generic machine
learning inference on heterogeneous treatment effects in randomized experiments, with an
application to immunization in India.* Tech. rep. National Bureau of Economic Research.

Wager, S. and S. Athey (2018). *Estimation and inference of heterogeneous treatment effects using random forests*. In: Journal of the American Statistical Association 113.523,
pp. 1228–1242.

This work was generously funded by the UK Medical Research Council (Grant #: MR/T04487X/1)

