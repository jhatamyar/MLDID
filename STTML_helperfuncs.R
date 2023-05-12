#'
#'Data cleaning and other helper functions 
#'Julia Hatamyar, 11 May 2022
#'
#'========================================
#'
#




did_datapreprocess <- function(data = data, Y = "tx_mi2"){
  
  ## rename variables 
  names(data)[names(data) == 'tyear'] <- 'G'
  names(data)[names(data) == 'year'] <- 'period'
  names(data)[names(data) == 'mun'] <- 'id'
  names(data)[names(data) == Y] <- 'Y' #should work 
  
  ### delete other outcome vars that are not Y 
  outcomes <- c("tx_mi2", "tx_ma2", "tx_mi_fet2", "tx_mmat2")
  outcomes <- outcomes[outcomes != Y]
  data <- data[,!colnames(data) %in% c(outcomes)]
  
  ## now need to drop mun with missing Y 
  data <- data %>% 
    group_by(id) %>%
    filter(!any(is.na(Y)))
  
  data <- data %>% 
    group_by(id) %>%
    filter(!any(G == 1996))
  
  #set those treated after the sample to zero (to act as control)
  data$G[data$G > 2008] <- 0
  #data <- data[data$period <]
  
  tlist <- unique(data$period)[order(unique(data$period))]
  
  # list of treated groups (by time) from smallest to largest
  glist <- unique(data$G, )[order(unique(data$G))]
  # Only the treated groups
  glist <- glist[glist>0]
  
  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)
  ## how many unique ID name (mun)
  n <- length(unique(data$id))
  
  # order dataset wrt idname and tname
  data <- data[order(data$id, data$period),]
  
  return(list(data=data, tlist = tlist, glist = glist, nT = nT, nG = nG, n = n))
}

# preprocess for the loop function
did_datapreprocess_sim <- function(data = dta){
  
  tlist <- unique(data$period)[order(unique(data$period))]
  
  # list of treated groups (by time) from smallest to largest
  glist <- unique(data$G, )[order(unique(data$G))]
  # Only the treated groups
  glist <- glist[glist>0]
  
  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)
  ## how many unique ID name (mun)
  n <- length(unique(data$id))
  
  # order dataset wrt idname and tname
  data <- data[order(data$id, data$period),]
  
  return(list(data=data, tlist = tlist, glist = glist, nT = nT, nG = nG, n = n))
}


### function for processing the output of the loop
### taken from https://github.com/bcallaway11/did/blob/master/R/process_attgt.R
process_attgt <- function(attgt.list) {
  nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
  
  # create vectors to hold the results
  group <- c()
  att <- c()
  att.se <- c()
  att.csa <- c()
  att.csa.se <- c()
  tt <- c()
  i <- 1
  
  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$TAU_hat
      att.se[i] <- attgt.list[[i]]$std.err.est
      att.csa[i] <- attgt.list[[i]]$attgt
      #att.csa.se[i] <- attgt.list[[i]]$attgt.se
      i <- i+1
    }
  }
  
  list(group=group, att=att, att.se=att.se, att.csa=att.csa, att.csa.se=att.csa.se, tt=tt)
}

### function for processing the output of the loop
### taken from https://github.com/bcallaway11/did/blob/master/R/process_attgt.R
process_attgt_sim <- function(attgt.list) {
  nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
  
  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1
  
  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      ## need to replace the pre-treatment years with zero 
      if ((attgt.list[[i]]$year - attgt.list[[i]]$group)<= -2 ){
        att[i] <- 0.00
      }
      else {
      att[i] <- attgt.list[[i]]$TAU.or
      }
      #att.csa.se[i] <- attgt.list[[i]]$attgt.se
      i <- i+1
    }
  }
  
  list(group=group, att=att, tt=tt)
}






### we aren't using this anymore, but it's here anyway 
run_interacted_regression <- function(affected, data, elist) {
  # Initialize a list to store the OLS summaries
  #OLS_summary_list <- vector("list", nperiods + 1)
  all_interaction_terms <- c()
  
  for (variable in affected) {
    
    
    # Initialize an empty vector to store the interaction terms
    interaction_terms <- c()
    
    # Loop through e.times and create the interaction terms, skipping the reference period (e.times = -1)
    for (i in elist) {
      if (i != -1) {
        if (i < 0) {
          interaction_term <- paste0(variable, ":e_", abs(i))
        } else if (i >= 0) {
          interaction_term <- paste0(variable, ":e", i)
        }
        interaction_terms <- c(interaction_terms, interaction_term)
      } 
    }
    # Add the interaction terms for the current variable to the list of all interaction terms
    all_interaction_terms <- c(all_interaction_terms, interaction_terms)
    
  }
  
  # Create the regression formula
  formula_text <- paste("dynamic.cate ~", paste(all_interaction_terms, collapse = " + "), "+ X1 + X2 + X3 + X4 + X5 + e_7 + e_6 +e_5 + e_4 + e_3 + e_2 + e0 + e1 + e2 + e3 + e4 + e5 + e6 + as.factor(period)")    
  # Run the regression
  OLSsim <- lm(formula_text, data = data)
  
  return(summary(OLSsim))
}

run_interacted_regressions <- function(nperiods, data) {
  # Initialize a list to store the OLS summaries
  OLS_summary_list <- vector("list", nperiods + 1)
  
  for (i in 0:nperiods) {
    # Create the interaction term
    interaction_term <- paste0("X1:e", i)
    
    # Create the regression formula
    formula_text <- paste("dynamic.cate ~", interaction_term, "+ X1 + e_3 + e_2 + e0 + e1 + e2")
    
    # Run the regression
    OLSsim <- lm(formula_text, data = data)
    
    # Store the summary in the list
    OLS_summary_list[[i + 1]] <- summary(OLSsim)
  }
  
  return(OLS_summary_list)
}



  

