### This code will be the master script for running simulations on Viking 
## Update from previous version to include CLANs and corrected BLP regressions
## as well as flexibility for different time periods 
### TO DO: get this to work with the "noise" extra 100 covariates version 


############################## PRELIMINARIES #########################

#Set to my WD
setwd("/Users/juliahatamyar/Documents/Research/BrazilML")
#setwd("~/scratch/STTML/SIMS29march")

#load required libraries
rm(list=ls(all=TRUE))
vec.pac= c("magrittr", "matrixStats", "Matrix", "BMisc", "grf", "did",
           "DRDID", "dplyr", "doParallel", "scales", "multcomp")

#install packages if not already installed 
list.of.packages <- vec.pac
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(vec.pac, require, character.only = TRUE)

#load functions 
source("STTML_helperfuncs.R")         ## pre and post processing 
source("STTML_ATTestimator.R")        ## the Lu, Nie, Wager estimator
source("STTML_aggregation.R")         ## Aggregation of group-times and CATEs/scores
source("STTML_heterogeneityfuncs.R")  ## Post-estimation (BLP, CLAN)
source("STTML_loop.R")                ## the main loop for Callaway Sant'Anna algorithm
source("generatesimdata.R")           ## data generating function 


### for saving the files,
## this changes according to naming convention in Google doc 
header <- "H1C0B0X0N2500."

## comment out if on Viking 
taskId <- (as.numeric("12"))

## uncomment if on Viking - This assigns the current array number to the taskID and seed
#taskIdChar <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#taskId <- (as.numeric(taskIdChar))

set.seed(round(taskId)) 


############## MAIN MLDID #####

## generate the data 
data.sim <- build_sim_datasetJH(time.periods = 4, n = 6000, 
                                random=TRUE, #change to false for confounding
                                noise = FALSE, #change to true to add 100 extra covariates
                                time.dependent.covar = TRUE, #adds X5
                                confounding = 1, #change to make confounding depend on more than 1 covar
                                het = "none", #this is heterogeneity in BETA NOT TAU
                                chi = 1, #how many covars have beta effect
                                taumodel = 1) 

#save the data (just in case as it has pscores)
write.csv(data, file = paste("data.", header, taskId, ".csv",sep=""))

#### run the main MLDID loop - this can take a long time 
res.sim <- didMLloop(data.sim,
                     xformla = ~X1 + X2 + X3 + X4 + X5)


## uncomment to save the results object 
##save(res.sim, file = paste("TESTRESULTS.", header, taskId, ".rda", sep="")) 

#### the main loop to get oracles 
res.sim.or <- didMLloop_oracle(data.sim,
                               xformla = ~X1 + X2 + X3 + X4 + X5)

## compare to vanilla callaway santanna, from the did R package
benchmark.sim <- att_gt(yname = "Y",
                        gname = "G",
                        idname = "id",
                        tname = "period",
                        xformla = ~X1 + X2 + X3 + X4 + X5,
                        data = data.sim,
                        #allow_unbalanced_panel=TRUE,
                        control_group = "notyettreated",
                        base_period = "universal",
                        est_method = "dr"
)


## extract the parameters for later 
parameters <- res.sim[["params"]]
## save number of groups and number of time periods 
n.groups <- parameters$nG
time.periods <- parameters$nT
# save the list of groups and periods  
tlist <- parameters$tlist
glist <- parameters$glist

## extract formatted data (I don't think I need to do this but just in case)
dta <- parameters$data

#### set up vector of times to event 
e.times <- sort(unique(data.sim$G - data.sim$period))



########### ATT(g,t) ################################################

# extract ATT(g,t) 
attgt.list.sim <- res.sim$attgt
attgt.list.CSA.sim <- benchmark.sim[["att"]]
attgt.list.or.sim <- res.sim.or$attgt

# process results
attgt.results.sim <- process_attgt(attgt.list.sim)
attgt.results.or.sim <- process_attgt_sim(attgt.list.or.sim)

#store for table 
sim.OR.attgt       <- attgt.results.or.sim$att
sim.MLDID.attgt    <- attgt.results.sim$att
sim.MLDID.attgt.se <- attgt.results.sim$att.se
sim.DRDID.attgt    <- benchmark.sim$att
sim.DRDID.attgt.se <- benchmark.sim$se

# Initialize an empty matrix for the table 
simtable.groupatt <- matrix(nrow = n.groups * time.periods, ncol = 3)

# Fill in with the appropriate values for each method 
for (i in 1:n.groups) {
  for (j in 1:time.periods) {
    row_index <- (i - 1) * time.periods + j
    simtable.groupatt[row_index, ] <- c(sim.OR.attgt[row_index],
                                        sim.MLDID.attgt[row_index],
                                        sim.DRDID.attgt[row_index])
  }
}


colnames(simtable.groupatt) <- c("OR", "MLDID", "DRDID")

# Initialize an empty matrix for the standarad errors 
simtable.groupatt.se <- matrix(nrow = n.groups * time.periods, ncol = 3)

# Fill in the matrix 
for (i in 1:n.groups) {
  for (j in 1:time.periods) {
    row_index <- (i - 1) * time.periods + j
    simtable.groupatt.se[row_index, ] <- c(0,
                                           sim.MLDID.attgt.se[row_index],
                                           sim.DRDID.attgt.se[row_index])
  }
}


## save the repetition as an R object 
save(simtable.groupatt, file = paste("groupatt.", header, taskId, ".rda", sep="")) 
save(simtable.groupatt.se, file = paste("groupatt.se.", header, taskId, ".rda", sep="")) 

##################### DYNAMIC ATTS #############################################

## aggregate the 3 versions 
ATT.dynamic.OR    <- aggregate_attgt_or(data=data.sim, attgt.res=attgt.results.or.sim)
ATT.dynamic.MLDID <- aggregate_attgt(data=data.sim, attgt.res=attgt.results.sim)
ATT.dynamic.DRDID <- aggte(benchmark.sim, type = "dynamic")

## create table 
simtable.dynamicATT <- cbind(ATT.dynamic.OR$att.e, ATT.dynamic.MLDID$dynamic.att.e, ATT.dynamic.DRDID$att.egt)
simtable.dynamicATT.se <- cbind(rep(0, length(ATT.dynamic.OR$att.e)), ATT.dynamic.MLDID$dynamic.att.e.se, ATT.dynamic.DRDID$se.egt)

colnames(simtable.dynamicATT) <- c("OR", "MLDID", "DRDID")
rownames(simtable.dynamicATT) <- as.character(e.times[2:(length(e.times) - 1)])

## save the run 
save(simtable.dynamicATT, file = paste("dynatt.", header, taskId, ".rda", sep="")) 
save(simtable.dynamicATT.se, file = paste("dynatt.se.", header, taskId, ".rda", sep="")) 


####################### DYNAMIC CATES AND SCORES #############################

## first, extract group-time CATEs
cates.gtt.OR    <- res.sim.or[["cates"]]
cates.gtt.MLDID <- res.sim[["cates"]]
scores.gtt.MLDID <- res.sim[["scores"]]
## saving raw cates just in case
save(cates.gtt.OR, file = paste("cates.gtt.OR.", header, taskId, ".rda", sep="")) 
save(cates.gtt.MLDID, file = paste("cates.gtt.MLDID.", header, taskId, ".rda", sep="")) 
save(scores.gtt.MLDID, file = paste("scores.gtt.MLDID.", header, taskId, ".rda", sep="")) 

## aggregate to dynamic (can be slow)
cates.dynamic.OR <- aggregate_cates(data=dta, res=res.sim.or, nG = n.groups, nT = time.periods)
cates.dynamic.OR <- as.data.frame(cates.dynamic.OR)

cates.dynamic.MLDID <- aggregate_cates(data=dta, res=res.sim, oracle = F, nG = n.groups, nT = time.periods)
cates.dynamic.MLDID <- as.data.frame(cates.dynamic.MLDID)

scores.dynamic.MLDID <- aggregate_scores(data=dta, res=res.sim, nG = n.groups, nT = time.periods)
scores.dynamic.MLDID <- as.data.frame(scores.dynamic.MLDID)

## saving them for RMSE later
save(cates.dynamic.OR, file = paste("cates.dyn.OR.", header, taskId, ".rda", sep="")) 
save(cates.dynamic.MLDID, file = paste("cates.dyn.MLDID.", header, taskId, ".rda", sep="")) 
save(scores.dynamic.MLDID, file = paste("scores.dyn.MLDID.", header, taskId, ".rda", sep="")) 


##################### BLP ANALYSIS ########################

## combine dynamic cates/scores with original data 
#### set up vector of time to events (modify as we estimate fewer than in list
e.times <- e.times[2:(length(e.times) - 1)]
dta$e.times <- dta$period - dta$G
## using the new function (very slow)
BLP.OR <- BLP_data(cates.dynamic.OR, dta, e.times)
## also remove untreated
BLP.OR <- BLP.OR[BLP.OR$G > 0,]
## true cate doesn't exist pre-treatment, so replace with zero 
BLP.OR$dynamic.cate <- ifelse(BLP.OR$e.times <0, 0, BLP.OR$dynamic.cate)

## Same for others (there should be a way to do this all at once but whatever)
BLP.cates <- BLP_data(cates.dynamic.MLDID, dta, e.times)
## also remove untreated
BLP.cates <- BLP.cates[BLP.cates$G > 0,]


BLP.scores <- BLP_data(scores.dynamic.MLDID, dta, e.times)
## also remove untreated
BLP.scores <- BLP.scores[BLP.scores$G > 0,]

## create list of variables to test with CLANS
affected <- c("X1", "X2", "X3", "X4", "X5")
affected_str <- paste(affected, collapse = " + ")

### run the pooled BLP (TO DO: make a function )
OLSsim.or <- lm(dynamic.cate ~ X1 + X2 + X3 + X4 + X5 + as.factor(period), data = BLP.OR)
#summary(OLSsim.or)
out.OLSsim.or <- summary(OLSsim.or)

OLSsim.cates <- lm(dynamic.cate ~ X1 + X2 + X3 + X4 + X5 + as.factor(period), data = BLP.cates)
#summary(OLSsim.cates)
out.OLSsim.cates <- summary(OLSsim.cates)

OLSsim.scores <- lm(dynamic.cate ~ X1 + X2 + X3 + X4 + X5 + as.factor(period), data = BLP.scores)
#summary(OLSsim.cates)
out.OLSsim.scores <- summary(OLSsim.scores)


# now by event times using the function : note nperiods should be 1 less than time periods
BLP.bye.OR <- BLP_eventtimes(data = BLP.OR, nperiods = 5, rhs_formula = affected_str)
BLP.bye.cates <- BLP_eventtimes(data = BLP.cates, nperiods = 5, rhs_formula = affected_str)
BLP.bye.scores <- BLP_eventtimes(data = BLP.scores, nperiods = 5, rhs_formula = affected_str)

### TO DO: save and run clans and save 
### 
save(BLP.bye.OR, file = paste("BLP.bye.OR.", header, taskId, ".rda", sep="")) 
save(BLP.bye.cates, file = paste("BLP.bye.cates.", header, taskId, ".rda", sep="")) 
save(BLP.bye.scores, file = paste("BLP.bye.scores.", header, taskId, ".rda", sep="")) 

#### CLANS

## first the glht version 
CLAN.OR     <- CLAN_glhtest(BLP.OR, affected, times = 4)
CLAN.cates  <- CLAN_glhtest(BLP.cates, affected, times = 4)
CLAN.scores <- CLAN_glhtest(BLP.scores, affected, times = 4)

## repeat for simple means test 
CLAN.OR.ttest     <- CLAN_ttest(BLP.OR, affected, time.periods = 4)
CLAN.cates.ttest  <- CLAN_ttest(BLP.cates, affected, time.periods = 4)
CLAN.scores.ttest <- CLAN_ttest(BLP.scores, affected, time.periods = 4)


## save both versions (we will report the glht version but just in case)
save(CLAN.OR, file = paste("CLAN.OR.", header, taskId, ".rda", sep="")) 
save(CLAN.cates, file = paste("CLAN.cates.", header, taskId, ".rda", sep="")) 
save(CLAN.scores, file = paste("CLAN.scores.", header, taskId, ".rda", sep="")) 

save(CLAN.OR.ttest, file = paste("CLAN.OR.ttest.", header, taskId, ".rda", sep="")) 
save(CLAN.cates.ttest, file = paste("CLAN.cates.ttest.", header, taskId, ".rda", sep="")) 
save(CLAN.scores.ttest, file = paste("CLAN.scores.ttest.", header, taskId, ".rda", sep="")) 


