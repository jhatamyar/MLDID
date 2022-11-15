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
  data$G[data$G > 2012] <- 0
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
  att.csa <- c()
  tt <- c()
  i <- 1
  
  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$TAU_hat
      att.csa[i] <- attgt.list[[i]]$attgt
      i <- i+1
    }
  }
  
  list(group=group, att=att, att.csa=att.csa, tt=tt)
}

### function for aggregating atts 
### takes output of above function as input 
aggregate_attgt <- function(data=data, attgt.res = attgt.res) {
## indicates year of first treatment 
  group <- attgt.res$group
## the group-time att 
  att <- attgt.res$att
## the group-time att (Callaway vversion)
  att.csa <- attgt.res$att.csa
## indicates current year 
  tt <- attgt.res$tt

# we can work in overall probabilities because conditioning will cancel out
# cause it shows up in numerator and denominator
  pg <- sapply(glist, function(g) mean(data$period==g))

# length of this is equal to number of groups
  pgg <- pg

# same but length is equal to the number of ATT(g,t)
## THIS WAS THE MISSING STEP OMG
  pg <- pg[match(group, glist)]

  originalt <- tt
  originalgroup <- group


## get the indicators for time to event (0 is event time)
  eseq <- unique(originalt - originalgroup)
  eseq <- eseq[order(eseq)]

# compute atts that are specific to each event time
  dynamic.att.e <- sapply(eseq, function(e) {
  # keep att(g,t) for the right g&t as well as ones that
  # are not trimmed out from balancing the sample
    whiche <- which( (originalt - originalgroup == e))
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })

# compute atts that are specific to each event time using CSA
  dynamic.att.e.csa <- sapply(eseq, function(e) {
  # keep att(g,t) for the right g&t as well as ones that
  # are not trimmed out from balancing the sample
    whiche <- which( (originalt - originalgroup == e))
    atte <- att.csa[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })
  list(dynamic.att.e = dynamic.att.e, dynamic.att.e.csa = dynamic.att.e.csa,
       originalt = originalt, originalgroup = originalgroup, pg = pg)
}