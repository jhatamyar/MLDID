#'
#' Aggregating ATTs and CATES/scores 
#'Julia Hatamyar, new file 11 May 2023, original code earlier version  
#'
#'========================================
#'
#

### function for aggregating atts 
### takes output of function "process_attgt" as input  
aggregate_attgt <- function(data=data, attgt.res = attgt.res) {
  ## indicates year of first treatment 
  group <- attgt.res$group
  ## the group-time att 
  att <- attgt.res$att
  ## the group-time att SEs 
  att.se <- attgt.res$att.se
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
  
  # compute standard errors that are specific to each event time (a bit shady)
  dynamic.att.e.se <- sapply(eseq, function(e) {
    # keep att(g,t) for the right g&t as well as ones that
    # are not trimmed out from balancing the sample
    whiche <- which( (originalt - originalgroup == e))
    atte.se <- att.se[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte.se*pge)
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
  list(dynamic.att.e = dynamic.att.e, dynamic.att.e.se = dynamic.att.e.se, dynamic.att.e.csa = dynamic.att.e.csa,
       originalt = originalt, originalgroup = originalgroup, pg = pg)
}

### function for aggregating atts 
### takes output of above function as input 
aggregate_attgt_or <- function(data=data, attgt.res = attgt.res) {
  ## indicates year of first treatment 
  group <- attgt.res$group
  ## the group-time att 
  att <- attgt.res$att
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
  dynamic.att.e.or <- sapply(eseq, function(e) {
    # keep att(g,t) for the right g&t as well as ones that
    # are not trimmed out from balancing the sample
    whiche <- which( (originalt - originalgroup == e))
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })
  list(att.e = dynamic.att.e.or, 
       originalt = originalt, originalgroup = originalgroup, pg = pg)
}

### function for aggregating CATES 
### takes output of above function as input 
aggregate_cates <- function(data=data, res = res, oracle = TRUE,
                            nG = nG, nT = nT) {
  
  ## first extract cates
  cates <- res$cates
  
  if (oracle == TRUE) {
    #cates <- cates[1:4020,]
    cates.agg <- do.call(rbind,
                         lapply(seq(1, nrow(cates), 2), function(i){
                           x <- cates[ i:(i + 1), , drop = FALSE]
                           cates.agg <- ifelse(is.na(x[1,]), x[2,], x[1,])
                         }))
    #cates.agg[cates.agg == 0] <- .
    
  }
  else {
    cates.agg <- do.call(rbind,
                         lapply(seq(1, nrow(cates), 2), function(i){
                           x <- cates[ i:(i + 1), , drop = FALSE]
                           cates.agg <- colSums(x)/2
                         }))
  }
  
  
  
  ## need oritinal groups and times to get event times out
  group <- c()
  tt <- c()
  i <- 1
  
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- res[["attgt"]][[i]][["group"]]
      tt[i] <- res[["attgt"]][[i]][["year"]]
      i <- i+1
    }
  }
  
  originalt <- tt
  originalgroup <- group
  pg <- sapply(glist, function(g) mean(data$period==g))
  
  pg <- pg[match(group, glist)]
  ## get the indicators for time to event (0 is event time)
  eseq <- unique(data$period - data$G)
  eseq <- eseq[order(eseq)]
  
  # create empty matrix
  dynamic.cates <- matrix(NA, nrow(cates.agg), length(eseq))
  
  for (i in 1:nrow(cates.agg)) {
    dynamic.cates[i,] <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (originalt - originalgroup == e))
      catte <- cates.agg[i,][whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(catte*pge, na.rm=TRUE)
      
    })
  }
  
  dynamic.cates <- cbind(unique(data$id), dynamic.cates)
  list(dynamic.cates = dynamic.cates)
}

### function for aggregating scores 
### takes output of above function as input 
aggregate_scores <- function(data=data, res = res, nG = nG, nT = nT) {
  
  ## first extract cates
  scores <- res$scores
  
  
  scores.agg <- do.call(rbind,
                        lapply(seq(1, nrow(scores), 2), function(i){
                          x <- scores[ i:(i + 1), , drop = FALSE]
                          scores.agg <- colSums(x)/2
                        }))
  
  ## need oritinal groups and times to get event times out
  group <- c()
  tt <- c()
  i <- 1
  
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- res[["attgt"]][[i]][["group"]]
      tt[i] <- res[["attgt"]][[i]][["year"]]
      i <- i+1
    }
  }
  
  originalt <- tt
  originalgroup <- group
  pg <- sapply(glist, function(g) mean(data$period==g))
  
  pg <- pg[match(group, glist)]
  ## get the indicators for time to event (0 is event time)
  eseq <- unique(data$period - data$G)
  eseq <- eseq[order(eseq)]
  
  # create empty matrix
  dynamic.scores <- matrix(NA, nrow(scores.agg), length(eseq))
  
  for (i in 1:nrow(scores.agg)) {
    dynamic.scores[i,] <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (originalt - originalgroup == e))
      scoree <- scores.agg[i,][whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(scoree*pge, na.rm = TRUE)
    })
  }
  dynamic.scores <- cbind(unique(data$id), dynamic.scores)
  list(dynamic.scores = dynamic.scores)
}  


