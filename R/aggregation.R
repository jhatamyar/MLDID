#'
#' Aggregating ATTs and CATES/scores 
#'Julia Hatamyar, new file 11 May 2023, original code earlier version  
#'
#'========================================
#'
#

#' @title Aggregate ATTs for Dynamic Event Study
#'
#' @description This function aggregates average treatment effects on the treated (ATTs) to create a dynamic 
#' "event study" type analysis. It takes an "MLDID" object as input, processes the results, and 
#' computes ATTs specific to each event time.
#'
#' @param attgt.result An object of class "MLDID" containing the ATT results to be aggregated.
#'
#' @return A list containing the processed "MLDID" object, dynamic ATT estimates, standard errors, 
#'         covariate-specific ATTs, and other relevant data used in the computation.
#'
#' @export
dynamic_attgt <- function(attgt.result) {
  ## process the results for clean format
  # Process the results and update the MLDID object
  attgt.result <- process_attgt(attgt.result)
  
  # Extract data and parameters from the MLDID object
  data <- attgt.result@params$data
  glist <- attgt.result@params$glist
  tlist <- attgt.result@params$tlist
  eseq <- attgt.result@params$e.times
  
  # Extract the processed results from the updated MLDID object
  processed_data <- attgt.result@processed_attgt
  group <- processed_data$group
  att <- processed_data$att
  att.se <- processed_data$att.se
  att.csa <- processed_data$att.csa
  tt <- processed_data$tt
  
  # we can work in overall probabilities because conditioning will cancel out
  # cause it shows up in numerator and denominator
  pg <- sapply(glist, function(g) mean(data$period==g))
  
  # length of this is equal to number of groups
  pgg <- pg
  
  # same but length is equal to the number of ATT(g,t)
  pg <- pg[match(group, glist)]
  
  originalt <- tt
  originalgroup <- group
  
  
  # compute atts that are specific to each event time
  dynamic.att.e <- sapply(eseq, function(e) {
    # keep att(g,t) for the right g&t as well as ones that
    # are not trimmed out from balancing the sample
    whiche <- which( (originalt - originalgroup == e))
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })
  
  # compute standard errors that are specific to each event time (a bit shady without IF)
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
  list(attgt.result = attgt.result, dynamic.att.e = dynamic.att.e, dynamic.att.e.se = dynamic.att.e.se, dynamic.att.e.csa = dynamic.att.e.csa,
       originalt = originalt, originalgroup = originalgroup, pg = pg)
}






#' @title Aggregate Oracle ATTs
#'
#' @description This function aggregates oracle average treatment effects on the treated (ATTs) based on 
#' the provided data. It is designed for processing and aggregating ATTs in a dynamic event study 
#' context, using oracle or known ATTs.
#'
#' @param data A data frame 
#' @param attgt.res An object containing oracle ATT results 
#'                  with group-time information.
#'
#' @return A list containing the aggregated oracle ATT estimates and other relevant data 
#'         used in the computation.
#'
#' @export
dynamic_attgt_or <- function(data=data, attgt.res = attgt.res) {
  ## indicates year of first treatment 
  group <- attgt.res$group
  ## the group-time att 
  att <- attgt.res$att
  ## indicates current year 
  tt <- attgt.res$tt
  ## need glist 
  glist <- sort(unique(data.or$G)[unique(data.or$G)>0])
  # we can work in overall probabilities because conditioning will cancel out
  # cause it shows up in numerator and denominator
  pg <- sapply(glist, function(g) mean(data$period==g))
  
  # length of this is equal to number of groups
  pgg <- pg
  
  # same but length is equal to the number of ATT(g,t)
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







#' @title Structure CATES Dynamically 
#'
#' @description This function structures Conditional Average Treatment Effects (CATES) dynamically based on 
#' time to events. It takes an object of class "MLDID" as input and processes it to compute dynamic CATES.
#'
#' @param res An object of class "MLDID" containing the results from which CATES will be structured dynamically.
#' @param type A character vector specifying the type of data to be structured. 
#'             Options are "cates" for Conditional Average Treatment Effects or "scores" for double-robust scores.
#'             Default is "cates".
#' @param oracle A logical flag indicating whether the function is being used for simulations 
#'               (oracle = TRUE) or real data analysis (oracle = FALSE). Default is FALSE.
#'
#' @return A data frame containing the dynamically structured CATES or scores, with each row 
#'         corresponding to a unique identifier from the input data and columns for each event time.
#'
#' @export
dynamic_cates <- function(res,   ## the output of the MLDID function 
                          type=c("cates","scores"), ## must specify 
                          oracle = FALSE) { ## set to true for simulations
  
  # Ensure input is an MLDID object
  if (!inherits(res, "MLDID")) {
    stop("Input must be an MLDID object.")
  }
  
  
  # Extract data and parameters from the MLDID object
  data <- res@params$data
  nG <- res@params$nG
  nT <- res@params$nT
  glist <- res@params$glist
  eseq <- res@params$e.times
  
  
  # Extract cates/scores based on the specified type
  cates <- if (type == "cates") res@cates else res@scores
  
  if (oracle == TRUE) {
    #cates <- cates[1:4020,]
    cates.agg <- do.call(rbind,
                         lapply(seq(1, nrow(cates), 2), function(i){
                           x <- cates[ i:(i + 1), , drop = FALSE]
                           cates.agg <- ifelse(is.na(x[1,]), x[2,], x[1,])
                         }))
    
  }
  else if (oracle == FALSE) {
    cates.agg <- do.call(rbind,
                         lapply(seq(1, nrow(cates), 2), function(i) {
                           x <- cates[i:min(i + 1, nrow(cates)), , drop = FALSE]
                           #print(x)  # Inspect the values of x
                           means <- colMeans(x, na.rm = TRUE)
                           #print(means)  # Inspect the calculated means
                           means
                         }))
  }
  
  
  
  ## need oritinal groups and times to get event times out
  group <- vector("numeric", length = nG * nT)
  tt <- vector("numeric", length = nG * nT)
  i <- 1
  
  
  for (f in 1:nG) {
    for (s in 1:nT) {
      if (!is.null(res@attgt[[i]])) {
        group[i] <- res@attgt[[i]]$group
        tt[i] <- res@attgt[[i]]$year
      }
      i <- i+1
    }
  }
  
  originalt <- tt
  originalgroup <- group
  pg <- sapply(glist, function(g) mean(data$period==g))
  
  pg <- pg[match(group, glist)]
  ## TESTING THIS DELETE IF NEEDED
  pg <- ifelse(tt >= group, pg, 0)
  
  weights = matrix(nrow = length(res@attgt), ncol = 1)
  for (i in 1:length(res@attgt)) {
    weights[i,] = (length(res@attgt[[i]]$taus)/2)/length(unique(data$id))
  }
  
  # create empty matrix
  dynamic.cates <- matrix(NA, nrow(cates.agg), length(eseq))
  
  for (i in 1:nrow(cates.agg)) {
    dynamic.cates[i,] <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (originalt - originalgroup == e))
      catte <- cates.agg[i,][whiche]
      #whichcate <- which((glist - unique_G[i]) >=0)
      #catte <- catte[whichcate]
      pge <- pg[whiche] /(sum(pg[whiche]))
      sum(catte*pge, na.rm=TRUE)
      
    })
  }
  
  dynamic.cates <- cbind(unique(data$id), dynamic.cates)
  return(as.data.frame(dynamic.cates))
}


