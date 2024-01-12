#'
#'Data cleaning and other helper functions
#'Julia Hatamyar, 11 May 2022
#'
#'========================================
#'
#

#' @description Function to process arguments passed to the main methods
#'
#' @import dplyr
#' @param outcome A vector of outcomes
#' @param group A vector indicating group membership
#' @param time A vector indicating time periods
#' @param id_name A vector indicating data ID (state, municipality, etc)
#' @param data A dataframe containing your data
#'
#' @return an intermediate processed step to be used within the MLDID function
#'
#' @export
did_datapreprocess <- function(outcome,
                               group,
                               time,
                               id_name,
                               data){


  ## rename variables
  names(data)[names(data) == group] <- 'G'
  names(data)[names(data) == time] <- 'period'
  names(data)[names(data) == id_name] <- 'id'
  names(data)[names(data) == outcome] <- 'Y' #should work


  ## now need to drop anything with missing Y
  ## First tell user how many are missing

  missing_Y <- sum(is.na(data$Y))
  if (missing_Y != 0) {
    warning(paste0("There are ", missing_Y, " missing Y values in the data"))


  original_ids <- length(unique(data$id))

  data <- data %>%
    group_by(id) %>%
    filter(!any(is.na(Y)))

  dropped_ids <- original_ids - length(unique(data$id))

  warning(paste0("Dropped ", dropped_ids, " groups with missing Y from the data"))

  }

  ## report any missings
  n_missing <- nrow(data) -sum(complete.cases(data))
  if (n_missing != 0) {
    warning(paste0("There are ", n_missing, " missing values in the data, consider correcting"))
  }

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

  #### set up vector of times to event
  e.times <- sort(unique(data$G - data$period))

  # order dataset wrt idname and tname
  data <- data[order(data$id, data$period),]

  return(list(data=data, tlist = tlist, glist = glist, nT = nT, nG = nG, n = n, e.times = e.times))
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


#' Process ATTGT Results
#'
#' This function processes the results stored in an "MLDID" object, extracting and organizing
#' Average Treatment Effects on the Treated (ATT) estimates and related statistics.Usually this
#' will be performed as an intermediate step within the aggregation function but it might be
#' useful to see the group-time ATTs on their own.
#'
#' @param attgt.result An object of class "MLDID" containing the results to be processed.
#'
#' @return An updated "MLDID" object that includes the processed ATT results and other
#'         relevant statistics, such as standard errors and covariate-specific ATTs.
#'
#' @examples
#' # Example usage:
#' # processed_result <- process_attgt(your_MLDID_object)
#'
#' @export
process_attgt <- function(attgt.result) {

  if (!inherits(attgt.result, "MLDID")) {
    stop("Input must be an MLDID object.")
  }

  # Extract parameters and attgt list
  nG <- attgt.result@params$nG
  nT <- attgt.result@params$nT
  attgt.list <- attgt.result@attgt

  #nG <- attgt.result$params$nG
  #nT <- attgt.result$params$nT

  ## get the atts
  #attgt.list <- attgt.result$attgt
  #nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  #nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))

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
  ## extract the attgts for input to aggregation function

  # Update or add new slots to the MLDID object
  # Assuming 'processed_data' is a new slot in your MLDID definition
  # This step may involve modifying the setClass definition to include new slots
  updated_result <- attgt.result
  updated_result@processed_attgt <- list(group = group, att = att, att.se = att.se, att.csa = att.csa, att.csa.se = att.csa.se, tt = tt)

  # Return the updated MLDID object
  return(updated_result)

  #list(group=group, att=att, att.se=att.se, att.csa=att.csa, att.csa.se=att.csa.se, tt=tt)
}




#' Process ATTGT Simulation Results
#'
#' This function is specifically tailored for processing simulation results stored in an "MLDID"
#' object, which has already been run through the non-simulated version. It focuses on extracting and organizing simulation-specific ATT estimates.
#'
#' @param attgt.result An updated object of class "MLDID" containing the simulation results to be processed.
#'
#' @return An updated "MLDID" object that includes the processed simulation-specific ATT results,
#'         along with the grouping and time period information.
#'
#' @examples
#' # Example usage:
#' # processed_sim_result <- process_attgt_sim(your_MLDID_object)
#'
#' @export
process_attgt_sim <- function(attgt.result) {

  # Extract parameters and attgt list
  nG <- attgt.result@params$nG
  nT <- attgt.result@params$nT
  attgt.list <- attgt.result@attgt

  #nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  #nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))

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

  updated_result <- attgt.result
  updated_result@processed_attgt <- list(group=group, att=att, tt=tt)
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
  # Generate e terms dynamically based on elist, skipping -1
  e_terms <- sapply(elist[elist != -1], function(i) {
    if (i < 0) {
      return(paste0("e_", abs(i)))
    } else {
      return(paste0("e", i))
    }
  }, simplify = TRUE)

  e_terms <- paste(e_terms, collapse = " + ")

  # Create the regression formula
  formula_text <- paste("dynamic.cate ~ ", paste(all_interaction_terms, collapse = " + "), "+ as.factor(period) + ", e_terms)
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





