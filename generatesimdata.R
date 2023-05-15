#'
#' Code for simulating data, adapted from CSA
#' 
#' Julia Hatamyar Began 31 August 2022
#' 
#' Additions Feb/March 2023: time varying X5
#' Additions APril 2023: noise covariates 
#' 
#' ==========================================
#' 

#setwd("/Users/juliahatamyar/Documents/Research/BrazilML")
#setwd("//userfs/jfh528/w2k/Research/BrazilML")

build_sim_datasetJH <- function(time.periods = 4, n = 2500, 
                                random=TRUE, #change to false for confounding
                                noise = FALSE, #change to true to add 100 extra covariates
                                time.dependent.covar = TRUE, #adds X5
                                confounding = 1, #change to make confounding depend on more than 1 covar
                                het = "none", #this is heterogeneity in BETA NOT TAU
                                chi = 1, #how many covars have beta effect
                                taumodel = 1) {
  
  ## set sim parameters using the function from did package
  ### TO DO: implement a version myself, as we don't use all the outputs
  sp <- did::reset.sim(time.periods = time.periods, n = n)
  
  #-----------------------------------------------------------------------------
  # build dataset
  #-----------------------------------------------------------------------------
  time.periods <- sp$time.periods
  nt <- sp$nt
  #bett <- sp_list$bett
  bett <- rep(0,time.periods) ## effect by X for treated potential outcome (default 0)
  thet=sp$thet
  nu <- sp$nu
  theu <- sp$theu
  betu <- rep(0,time.periods) ## effect by X for untreated (default 0)
  te.bet.ind <- sp$te.bet.ind ## for now, = a vector of 1s, overall TE is 1
  te.bet.X <- bett # for now, =  bett 
  te.t <- sp$te.t ## for now, te.t = thet (no time effects)
  te.e <- sp$te.e ## default is zero, no dynamic effects 
  te <- sp$te
  n <- n
  gamG <- sp$gamG ## propensity score param: in the reset function this is going to be 0.5*g/G
  #reg <- sp_list$reg
  
  te.e <- 1:time.periods
  #ipw <- ipw 
  
  # generate covars
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)
  X4 <- rbinom(n, 1, 0.5) ## will not be used for anything, just noise
  
  ### generate time varying covar 
  if (time.dependent.covar) {
    # Generate a matrix of time points, with each row representing one observation
    time_points <- matrix(rep(seq(0, 1, length.out = time.periods), n), nrow = n, byrow = TRUE)  
    # Define a function that maps time points to covariate values
    covariate_function <- function(time) {
      sin(2 * pi * time) + rnorm(n)
    }
    
    # Evaluate the covariate function at each time point for each observation
    covariate <- apply(time_points, 2, covariate_function)
    
    # Convert the covariate matrix to a data frame with four columns for each observation
    X5tdf <- data.frame(matrix(covariate, nrow = n))
    
    # Name the columns of the time-dependent covariate data frame
    colnames(X5tdf ) <- paste0("X5_", 1:time.periods)
  } else { 
    X5tdf <- NULL
  }
  # now add the optional 100 "noisy covaraites
  if (noise) {
    noise_continuous <- matrix(rnorm(n * 50), nrow = n) # 50 continuous noise covariates
    noise_binary <- matrix(rbinom(n * 50, 1, 0.5), nrow = n) # 50 binary noise covariates
    noise_covariates <- cbind(noise_continuous, noise_binary) # Combine continuous and binary noise covariates
    colnames(noise_covariates) <- paste0("noise_", 1:100) # Set column names for noise covariates
  } else {
    noise_covariates <- NULL
  }
  
  ## set confounding complexity
  if (confounding == 1){
    c.model <- X2
  } else {
    c.model <- X1 + X2 + X3 
  }
  
  ## set heterogeneity params 
  if (het == "none"){
    te.e <- 1:time.periods
  }
  if (het == "constant"){
    te.e <- 1:time.periods
    te.bet.X <- rep(1, time.periods)
  }
  if (het == "calendar"){
    te.e <- 1:time.periods
    te.bet.X <- 1:time.periods
  }
  
  ## print values for my sanity
  cat(paste("current bett:", bett, "\n"))
  cat(paste("current betu:", betu, "\n"))
  cat(paste("current te.e:", te.e, "\n"))
  cat(paste("current te.bet.X:", te.bet.X, "\n"))
  #cat(paste("current ipw:", ipw, "\n"))
  cat(paste("current het:", het, "\n"))
  cat(paste("current te:", te, "\n"))
  
  # making pr depend on covariates 
  if (random == TRUE) {
    pr <- matrix(1 / (time.periods + 1), nrow = n, ncol = (time.periods + 1)) ## this is random assignment 
  } else {
    pr <- exp(outer(c.model ,gamG)) / apply( exp(outer(c.model,gamG)), 1, sum)
  }
  
  G <- apply(pr, 1, function(pvec) sample(seq(0,time.periods), size=1, prob=pvec))
  
  Gt <- G[G>0] ## gets the treated units 
  nt <- length(Gt)
  
  ## updating this to include all covars
  if (chi == 1) {
    Xmodel <- X1 
  } else {
    Xmodel <- (X1 + X2 + X3)
  }
  
  Xt <- Xmodel[G>0] # X for the treated group 
  
  # draw individual fixed effect
  Ct <- rnorm(nt, mean=G)
  
  # generate untreated potential outcomes in each time period
  Ynames <- paste0("Y",1:time.periods)
  
  Y0tmat <- sapply(1:time.periods, function(t) {
    thet[t] + Ct + Xt*bett[t] + rnorm(nt)
  })
  Y0tdf <- as.data.frame(Y0tmat)
  
  ## names for untreated potential outcomes
  Y0names <- paste0("y0_",1:time.periods)
  colnames(Y0tdf) <- Y0names
  
  # generate treated potential outcomes
  
  ## NOW MAKING THE TAUS 
  if (taumodel == 1) {
    tau <- X1[G>0]
  } else {
    tau <- (X2[G>0] + X3[G>0])^2
  }
  
  cat(paste("The value of random within the function is:", random, "\n"))
  
  Y1tdf <- sapply(1:time.periods, function(t) {
    te.t[t] + te.bet.ind[Gt]*Ct + Xt*te.bet.X[t] 
    + ((Gt <= t)*te.e[sapply(1:nt, function(i) max(t-Gt[i]+1,1))] ## default of this is zero, but we added dynamic effects 
       + te)*tau + rnorm(nt) #
  })
  
  ## names for treated potential outcomes
  Y1names <- paste0("y1_",1:time.periods)
  colnames(Y1tdf) <- Y1names
  
  # generate observed data
  Ytdf <- sapply(1:time.periods, function(t) {
    (Gt<=t)*Y1tdf[,t] + (Gt>t)*Y0tdf[,t]
  })
  colnames(Ytdf) <- Ynames
  
  ## combine unobserved data 
  # store observed data for treated group
  #dft <- cbind.data.frame(G=Gt,X1=X1[G>0],X2=X2[G>0],X3=X3[G>0],X4=X4[G>0],X5tdf[G>0,],Ytdf,Y0tdf,Y1tdf)
  #dft <- cbind.data.frame(G = Gt, X1 = X1[G > 0], X2 = X2[G > 0], X3 = X3[G > 0], X4 = X4[G > 0], 
  # X5tdf[G > 0, ], Ytdf, Y0tdf, Y1tdf, noise_covariates[G > 0, ])
  
  ## combine unobserved data
  if (time.dependent.covar & noise) {
    dft <- cbind.data.frame(G = Gt, X1 = X1[G > 0], X2 = X2[G > 0], X3 = X3[G > 0], X4 = X4[G > 0], X5tdf[G > 0, ], Ytdf, Y0tdf, Y1tdf, noise_covariates[G > 0, ])
  } else if (time.dependent.covar & !noise) {
    dft <- cbind.data.frame(G = Gt, X1 = X1[G > 0], X2 = X2[G > 0], X3 = X3[G > 0], X4 = X4[G > 0], X5tdf[G > 0, ], Ytdf, Y0tdf, Y1tdf)
  } else if (!time.dependent.covar & noise) {
    dft <- cbind.data.frame(G = Gt, X1 = X1[G > 0], X2 = X2[G > 0], X3 = X3[G > 0], X4 = X4[G > 0], Ytdf, Y0tdf, Y1tdf, noise_covariates[G > 0, ])
  } else {
    dft <- cbind.data.frame(G = Gt, X1 = X1[G > 0], X2 = X2[G > 0], X3 = X3[G > 0], X4 = X4[G > 0], Ytdf, Y0tdf, Y1tdf)
  }
  
  
  ## also need probability of being treated in times after t=1
  Gprobs <- pr[,3:(length(seq(0,time.periods)))]
  ## setting colnames conditionally so we can have a version with only 2 periods
  if (time.periods == 2){
    Gprobs <- matrix(Gprobs, ncol = 1)
    colnames(Gprobs) <- c("G2")
  } else {
    colnames(Gprobs) <- paste0("G", 2:time.periods)
  }
  # untreated units
  
  # draw untreated covariate
  nu <- sum(G==0)
  Xu <- Xmodel[G==0]
  
  # draw untreated fixed effect
  Cu <- rnorm(nu, mean=0)
  
  
  # generate untreated potential outcomes
  Y0umat <- sapply(1:time.periods, function(t) {
    theu[t] + Cu + rnorm(nu) + Xu*betu[t]
  })
  Y0udf <- as.data.frame(Y0umat)
  colnames(Y0udf) <- Ynames
  
  
  y0udf <- Y0udf
  colnames(y0udf) <- Y0names
  y1udf <- Y0udf
  colnames(y1udf) <- Y1names
  
  # store dataset of observed outcomes for untreated units
  #dfu <- cbind.data.frame(G=0,X1=X1[G==0],X2=X2[G==0],X3=X3[G==0],X4=X4[G==0],X5tdf[G==0,], Y0udf, y0udf,y1udf)
  #dfu <- cbind.data.frame(G = 0, X1 = X1[G == 0], X2 = X2[G == 0], X3 = X3[G == 0], X4 = X4[G == 0], 
  # X5tdf[G == 0, ], Y0udf, y0udf, y1udf, noise_covariates[G == 0, ])
  
  if (time.dependent.covar & noise) {
    dfu <- cbind.data.frame(G = 0, X1 = X1[G == 0], X2 = X2[G == 0], X3 = X3[G == 0], X4 = X4[G == 0], X5tdf[G == 0, ], Y0udf, y0udf, y1udf, noise_covariates[G == 0, ])
  } else if (time.dependent.covar & !noise) {
    dfu <- cbind.data.frame(G = 0, X1 = X1[G == 0], X2 = X2[G == 0], X3 = X3[G == 0], X4 = X4[G == 0], X5tdf[G == 0, ], Y0udf, y0udf, y1udf)
  } else if (!time.dependent.covar & noise) {
    dfu <- cbind.data.frame(G = 0, X1 = X1[G == 0], X2 = X2[G == 0], X3 = X3[G == 0], X4 = X4[G == 0], Y0udf, y0udf, y1udf, noise_covariates[G == 0, ])
  } else {
    dfu <- cbind.data.frame(G = 0, X1 = X1[G == 0], X2 = X2[G == 0], X3 = X3[G == 0], X4 = X4[G == 0], Y0udf, y0udf, y1udf)
  }
  
  
  # store overall datasetof observed data
  df <- rbind.data.frame(dft, dfu)
  
  # generate id variable
  df$id <- 1:nrow(df)
  # generate clusters (there's no actual within-cluster correlation)
  df$cluster <- sample(1:50, size=nrow(df), replace=TRUE)
  
  #cat(paste("current te.e:", te.e, "\n"))
  
  # convert data from wide to long format
  ddf <- tidyr::pivot_longer(df,
                             cols=c(tidyr::starts_with("X5_"),tidyr::starts_with("Y"),tidyr::starts_with("y0_"),tidyr::starts_with("y1_")),
                             names_to=c(".value", "period"),
                             names_pattern="(.*)(.)",
                             #names_prefix=c("Y", "y0", "y1"),
                             #names_sep = 1,
                             values_to=c("Yobs", "y0", "y1"))
  
  colnames(ddf)[colnames(ddf) == "X5_"] <- "X5"
  
  
  ddf$period <- as.numeric(ddf$period)
  ddf$treat <- 1*(ddf$G > 0)
  ## also need event time 
  ddf$delta.e <- ddf$period - ddf$G + 1 
  ddf$delta.e <- ifelse(ddf$delta.e < 0, 0, ddf$delta.e)
  ddf <- cbind(ddf, Gprobs)
  ddf <- ddf[order(ddf$id, ddf$period),] # reorder data
  
  ddf <- subset(ddf, G != 1)
  ddf
}

