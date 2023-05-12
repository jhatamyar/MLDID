### This file contains the helper functions for post-estimation inference
### AKA the treatment effect heterogeneity functions for BLP and CLAN



## function that merges dynamic cates/scores to the original dataset 
### by event time 
BLP_data <- function(cates=cates, data = data, e.times = e.times) {
  cate.df <- list()
  ### loop through dynamic cates 
  for (i in 1:nrow(cates)) {
    ## get one row 
    cate.temp <- as.data.frame(cates[i,])
    cate.temp.vals <- as.matrix(cate.temp[,2:(length(e.times)+1)])
    ## get the mun ID 
    dismun <- cate.temp[1]
    ## add event times 
    discate <- as.data.frame(cbind(t(cate.temp.vals), e.times))
    discate$id <- as.integer(dismun)
    colnames(discate) <- c("dynamic.cate", "e.times", "id")
    rownames(discate) <- NULL
    ## separate data, this is too slow
    datag <- data[data$id == as.integer(dismun),]
    ## merge together and save 
    temp.data <- merge(datag,discate, by = "e.times", all.x = TRUE)
    ### get the cate column and append 
    cate.df <- rbind(cate.df, temp.data)
    
  }
  # Loop through the integers and create new columns
  for (i in e.times) {
    column_name <- paste0("e", ifelse(i < 0, "_", ""), abs(i))
    cate.df[[column_name]] <- ifelse(cate.df$e.times == i, 1, 0)
  }
  
  return(cate.df)
}


## alternative to the above, need to test to see which is faster 
BLP_datav2 <- function(cates = cates, data = data, e.times = e.times) {
  process_row <- function(cate_row) {
    cate.temp <- as.data.frame(cate_row)
    cate.temp.vals <- as.matrix(cate.temp[, 2:(length(e.times) + 1)])
    
    dismun <- cate.temp[1]
    discate <- as.data.frame(cbind(t(cate.temp.vals), e.times))
    discate$id <- as.integer(dismun)
    colnames(discate) <- c("dynamic.cate", "e.times", "id")
    rownames(discate) <- NULL
    
    datag <- data[data$id == as.integer(dismun), ]
    temp.data <- merge(datag, discate, by = "e.times", all.x = TRUE)
    return(temp.data)
  }
  
  cate.df <- do.call(rbind, lapply(1:nrow(cates), function(i) process_row(cates[i, ])))
  return(cate.df)
}


## older version, only works for simulations 
BLP_eventtimes_old <- function(data = data, nperiods = 2) {
  
  # Initialize a list to store the OLS summaries
  OLSetimes_summary_list <- vector("list", (nperiods))
  
  for (i in 0:(nperiods-1)) {
    # Create temporary data frame for each e.times value
    data.temp <- data[data$e.times == i,]
    
    if (i < nperiods-1) {
      # Fit the linear model for each e.times value
      OLSsim <- lm(dynamic.cate ~ 1+  X1 + X2 + X3 + X4 + X5 + as.factor(period), data = data.temp)
    } else {
      ## must do this because the last event time only has one period 
      OLSsim <- lm(dynamic.cate ~ 1 + X1 + X2 + X3 + X4 + X5, data = data.temp)  
    }
    # Store the summary of the linear model in the list
    OLSetimes_summary_list[[i + 1]] <- summary(OLSsim)
    
  }  
  return(OLSetimes_summary_list)
  
  
}


### version with flexible RHS formula (so user can input whatever covaraites they want)
## note the formula must be a string object 
BLP_eventtimes <- function(data = data, nperiods = 2, rhs_formula = "X1 + X2 + X3 + X4 + X5") {
  
  # Initialize a list to store the OLS summaries
  OLSetimes_summary_list <- vector("list", (nperiods))
  
  for (i in 0:(nperiods-1)) {
    # Create temporary data frame for each e.times value
    data.temp <- data[data$e.times == i,]
    
    # Add "dynamic.cate ~ " to the beginning of the formula string
    lhs_formula <- "dynamic.cate ~ 1+ "
    full_formula <- paste(lhs_formula, rhs_formula)
    
    if (i < nperiods-1) {
      # Fit the linear model for each e.times value
      # Add as.factor(period) to the formula string
      full_formula_period <- paste(full_formula, "+ as.factor(period)")
      OLSsim <- lm(as.formula(full_formula_period), data = data.temp)
    } else {
      ## must do this because the last event time only has one period 
      OLSsim <- lm(as.formula(full_formula), data = data.temp)  
    }
    # Store the summary of the linear model in the list
    OLSetimes_summary_list[[i + 1]] <- summary(OLSsim)
    
  }  
  return(OLSetimes_summary_list)
  
}

##testing
#test <- lm(dynamic.cate ~ 1+ X1 + X2 + X3 + X4 + X5 + as.factor(period), data = data.temp.cates)
#summary(test)
### Adapted from Demirer code, as in Chernozhukov 2018 
CLAN_glhtest <- function(data = data, affected = affected, times = 2, thres = 0.2, alpha = 0.05){
  out.array <- array(NA, dim = c(length(affected), 15, times+1))
  
  for (i in 0:times) {
    # Create temporary data frame for each e.times value
    data.temp <- data[data$e.times == i,]
    
    # get the quantile based on the threshold 
    high.effect     <- quantile(data.temp$dynamic.cate, 1-thres)
    low.effect      <- quantile(data.temp$dynamic.cate, thres)
    
    # create indicator for being in top or bototm 
    data.temp$h       <- as.numeric(data.temp$dynamic.cate>high.effect)
    data.temp$l       <- as.numeric(data.temp$dynamic.cate<low.effect)
    
    for (m in 1:length(affected)) {
      
      form  <- paste(affected[m],"~h+l-1", sep="")
      reg   <- lm(form, data=data.temp[(data.temp$h==1)| (data.temp$l==1),])
      coef  <- reg$coefficients['h'] - reg$coefficients['l']
      test  <- glht(reg, linfct = c("h-l==0"))
      
      coef  <- (summary(reg)$coefficients['h',1])
      pval  <- (summary(reg)$coefficients['h',4])
      ## contains: coefficient, high/low CIs, then one sided test pvals 
      res1  <- c(coef, confint(reg, 'h', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      
      coef  <- (summary(reg)$coefficients['l',1])
      pval  <- (summary(reg)$coefficients['l',4])
      res2  <- c(coef, confint(reg, 'l', level = 1-alpha)[1:2], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      
      coef  <- (summary(reg)$coefficients['h',1]) - (summary(reg)$coefficients['l',1])
      pval  <- summary(test)$test$pvalues[1]
      res3  <- c((confint(test,level = 1-alpha))$confint[1:3], (as.numeric(coef<0)*(pval/2) + as.numeric(coef>0)*(1-pval/2)),(as.numeric(coef<0)*(1-pval/2) + as.numeric(coef>0)*(pval/2)))
      a     <- c(res1, res2, res3)
      
      out.array[m,,i+1] <- a
    }
    
  } 
  
  return(out.array)
  
}

CLAN_ttest <- function(data = data, affected = affected, time.periods = 3) {
  ## contains: coefficients, high/low CIs, ttest pval
  ## each row should be a different X
  out.array <- array(NA, dim = c(length(affected), 6, time.periods))
  thres <- 0.2
  ## loop through event times 
  for (i in 0:(time.periods - 1)) {
    # Create temporary data frame for each e.times value
    data.temp <- data[data$e.times == i,]
    
    # get the quantile based on the threshold 
    high.effect     <- quantile(data.temp$dynamic.cate, 1-thres)
    low.effect      <- quantile(data.temp$dynamic.cate, thres)
    
    # create indicator for being in top or bototm 
    data.temp$h       <- as.numeric(data.temp$dynamic.cate>high.effect)
    data.temp$l       <- as.numeric(data.temp$dynamic.cate<low.effect)
    
    # Split the dataframe into high and low effect groups
    high_effect_data <- data.temp[data.temp$h == 1,]
    low_effect_data <- data.temp[data.temp$l == 1,]
    
    # Initialize empty vectors to store results
    high_means <- numeric(length(affected))
    low_means <- numeric(length(affected))
    p_values <- numeric(length(affected))
    conf_ints <- matrix(NA, length(affected), 2)
    
    # Loop through the affected variables
    for (m in seq_along(affected)) {
      variable <- affected[m]
      
      # Calculate the means for high and low effect groups
      high_means[m] <- mean(high_effect_data[[variable]], na.rm = TRUE)
      low_means[m] <- mean(low_effect_data[[variable]], na.rm = TRUE)
      
      # Perform a t-test to check for significant differences
      test_result <- t.test(high_effect_data[[variable]], low_effect_data[[variable]], na.rm = TRUE)
      p_values[m] <- test_result$p.value
      conf_ints[m,] <- test_result$conf.int
      
      ## combine to matrix
      testres <- cbind(high_means, low_means, (high_means - low_means), conf_ints, p_values)
    }
    ## store to array 
    out.array[1:m,,i+1] <- testres
  }
  return(out.array)
} 

