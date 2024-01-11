### This file contains the helper functions for post-estimation inference
### AKA the treatment effect heterogeneity functions for BLP and CLAN



#' @title Prepare Heterogeneous Treatment Effects for Analysis
#'
#' @description Merges dynamic conditional average treatment effects (CATEs) or scores
#' with the original dataset by event-time. It is  used for conducting BLP
#' or CLAN analyses.
#'
#' @param attgt.result The output of the MLDID function.
#' @param cates The output of the dynamic_cates function, which could be CATEs or scores.
#'
#' @return A panel data frame where each row represents a unit (e.g., an individual or a firm)
#' at a specific event time, merged with the corresponding CATE or score.
#'
#' @export
het_prep <- function(attgt.result,  ## the output of MLDID
                     cates) {  ## the output of dynamic_cates (could be scores)
  data <- attgt.result@params$data
  e.times <- attgt.result@params$e.times

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
    ## need to create this in the data to merge on
    datag$e.times <- datag$period - datag$G
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


#' BLP Event Times Analysis
#'
#' This function facilitates BLP analysis with a flexible
#' RHS formula, allowing users to input any desired covariates. The function iteratively
#' fits a linear model for different event times (`e.times`) in the data.
#'
#' The formula for the linear model includes an intercept (denoted by `1` in the formula). Users must
#' provide the RHS formula as a string, which will be appended to `dynamic.cate ~ 1 +` to
#' form the full model formula.
#'
#' @param data The dataset to be used in the analysis, output from the function het_prep
#' @param nperiods The number of periods to consider in the analysis. Note this must be fewer than the number
#' of periods where e >=0
#' @param rhs_formula A string representing the RHS of the linear model formula, including
#' desired covariates. Default is NULL.
#'
#' @return A list of `nperiods` elements, each containing the summary of the linear model
#' fitted for a specific event time.
#'
#' @examples
#' # Assuming data is predefined and has appropriate structure:
#' results_list <- BLP_eventtimes(data, nperiods = 3, rhs_formula = "X1 + X2 + X3")
#'
#' @export
BLP_eventtimes <- function(data = data, nperiods = 2, rhs_formula = NULL) {

  # Initialize a list to store the OLS summaries
  OLSetimes_summary_list <- vector("list", (nperiods))

  for (i in 0:(nperiods)) {
    # Create temporary data frame for each e.times value
    data.temp <- data[data$e.times == i,]

    # Determine the formula based on rhs_formula
    if (is.null(rhs_formula)) {
      # If rhs_formula is NULL, use only the intercept
      full_formula <- "dynamic.cate ~ 1"
    } else {
      # Add "dynamic.cate ~ 1 + " to the beginning of the formula string
      lhs_formula <- "dynamic.cate ~ 1 + "
      full_formula <- paste(lhs_formula, rhs_formula)
    }

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


#' @title Generalized Linear Hypothesis Test in CLAN Analysis
#'
#' @description Adapted from the methodology in Demirer and Chernozhukov (2018), this function
#' performs a generalized linear hypothesis (GLH) test for each event time in the CLAN
#' (Classification Analysis) framework. It computes effects, confidence intervals,
#' and one-sided test p-values for specified variables based on a quantile-based
#' threshold, using quantiles of the estimated CATT.
#'
#' @import stats
#' @param data The dataset containing variables for analysis, output of function het_prep
#' @param affected A vector of variables for which the effects are to be tested.
#' @param thres A threshold value used to determine the high and low effect based on
#' dynamic.cate quantiles.
#' @param alpha Significance level used for confidence intervals.
#'
#' @return A three-dimensional array with dimensions corresponding to the number of
#' affected variables, the number of statistical measures (coefficients, CIs, one-sided p-values),
#' and the number of event times.
#'
#' @examples
#' # Assuming data and affected are predefined:
#' results_array <- CLAN_glhtest(data, affected, thres = 0.2, alpha = 0.05)
#'
#' @export
CLAN_glhtest <- function(data = data, affected = affected, thres = 0.2, alpha = 0.05){

  data <- data[data$G !=0,]

  times <- max(data$e.times)
  # Define the column names for the output array
  column_names <- c("coef_h", "CIlower_h", "CIupper_h", "pval_high", "pval_low",
                    "coef_l", "CIlower_l", "CIupper_l", "pval_low", "pval_high",
                    "coef_h-l", "CIlower_h-l", "CIupper_h-l", "pval_h-l_low", "pval_h-l_high")

  out.array <- array(NA, dim = c(length(affected), length(column_names), times+1), dimnames = list(affected, column_names, NULL))

 # out.array <- array(NA, dim = c(length(affected), 15, times+1))

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
      # Debug: Print current affected variable
      #print(paste("Processing affected variable:", affected[m]))

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

      # Debug: Check array assignment
      #print(paste("Storing results for variable", affected[m], "at time", i))


      out.array[m,,i+1] <- a
    }

  }

  return(out.array)

}


#' @title T-Test Analysis for CLAN
#'
#' @description This function, part of the CLAN framework, performs
#' t-tests to compare the means of affected variables between most and least effected groups.
#' These groups are determined based on dynamic.cate quantiles. The function operates
#' for each event time separately.
#'
#' @import stats
#'
#' @param data The dataset containing variables for analysis.
#' @param affected A vector of variables for which the t-tests are to be conducted.
#'
#' @return A three-dimensional array with dimensions corresponding to the number of
#' affected variables, the number of statistical measures (means, CIs, p-values),
#' and the number of event times.
#'
#' @export
CLAN_ttest <- function(data = data, affected = affected) {
  ## contains: coefficients, high/low CIs, ttest pval
  data <- data[data$G !=0,]
  ## each row should be a different X
  time.periods <- max(data$e.times) + 1

  column_names <- c("mean_h", "mean_l", "h-l", "CI_upper", "CI_lower", "pval")


  out.array <- array(NA, dim = c(length(affected), 6, time.periods), dimnames = list(affected, column_names, NULL))
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

#' Summary of BLP Model Outputs
#'
#' This function generates a summary table for the output of BLP_eventtimes
#' model analyses. It compiles the coefficients and their standard errors for each variable
#' across different event times into a single matrix.
#'
#' @param BLP_output A list containing the output of BLP_eventtimes
#' @param affected A character vector of the names of the variables included in the BLP model
#' analysis.
#'
#' @return A matrix where each row corresponds to an event time and each column represents
#' the coefficient or standard error of a variable. The matrix includes columns for both
#' coefficients and standard errors for each variable including the intercept.
#'
#' @examples
#' # Assuming BLP_output and affected are predefined:
#' summary_matrix <- BLPsummary(BLP_output, affected)
#' print(summary_matrix)
#'
#' @export
BLPsummary <- function(BLP_output, affected) {
  # Number of event times
  ntimes <- length(BLP_output)

  # Prepare column names for coefficients and standard errors
  variables <- c("(Intercept)", affected)
  colnames_coeff <- paste0(variables, "_coef")
  colnames_se <- paste0(variables, "_se")
  colnames_table <- c(colnames_coeff, colnames_se)

  # Initialize a matrix to store the coefficients and standard errors
  results_matrix <- matrix(NA, nrow = ntimes, ncol = length(colnames_table))
  colnames(results_matrix) <- colnames_table
  rownames(results_matrix) <- paste0("EventTime", 0:(ntimes-1))

  # Extract coefficients and standard errors for each event time
  for (i in 1:ntimes) {
    coeffs <- BLP_output[[i]]$coefficients[, 1]  # Extract coefficients
    ses <- BLP_output[[i]]$coefficients[, 2]     # Extract standard errors

    # Ensure the order of coefficients matches the order of 'variables'
    coeffs_ordered <- coeffs[match(variables, names(coeffs))]
    ses_ordered <- ses[match(variables, names(ses))]

    # Fill the matrix with coefficients and standard errors
    results_matrix[i, colnames_coeff] <- coeffs_ordered
    results_matrix[i, colnames_se] <- ses_ordered
  }

  return(results_matrix)
}




#' @title Plot BLP Coefficients Over Event Time
#'
#' @description Creates a plot of coefficient estimates and their confidence intervals
#' for various variables over different event times. The plot is generated based on
#' the output from the BLP_eventtimes function. It excludes the
#' intercept and presents each variable with a slight horizontal offset for readability.
#'
#' @import ggplot2
#'
#' @param BLP_output A list containing the output of `BLP_eventtimes`.
#' @param affected A character vector of the names of the variables included in the BLP model
#' analysis (excluding the intercept).
#'
#' @return A ggplot object representing the coefficient estimates and their confidence intervals
#' for each variable across different event times. Points are slightly dodged for readability,
#' and the plot displays integer values for event times on the x-axis.
#'
#' @examples
#' # Assuming BLP_output and affected are predefined:
#' plot <- plotBLP(BLP_output, affected)
#' plot
#'
#' @export
plotBLP <- function(BLP_output, affected) {
  plot_data <- data.frame(EventTime = integer(),
                          Estimate = numeric(),
                          LowerCI = numeric(),
                          UpperCI = numeric(),
                          Variable = factor())

  # Loop through each BLP output
  for (i in 1:length(BLP_output)) {
    coeffs_df <- BLP_output[[i]]$coefficients
    vars <- rownames(coeffs_df)

    # Debug: Check for NAs in vars
    #print(paste("Variables at Event Time", i, ":", paste(vars, collapse = ", ")))


    # Calculate confidence intervals
    # Assuming a 95% confidence level, z-value is approximately 1.96 for a normal distribution
    z_val <- 1.96
    lower_ci <- coeffs_df[, "Estimate"] - z_val * coeffs_df[, "Std. Error"]
    upper_ci <- coeffs_df[, "Estimate"] + z_val * coeffs_df[, "Std. Error"]

    # Combine data for plot
    for (var in vars) {
      if (var != "(Intercept)" && !is.na(var)){  # Exclude the intercept
        plot_data <- rbind(plot_data, data.frame(
          EventTime = i - 1,
          Estimate = coeffs_df[var, "Estimate"],
          LowerCI = lower_ci[var],
          UpperCI = upper_ci[var],
          Variable = var
      ))
    }
    }
  }

  # Drop unused factor levels and Filter out rows with NAs
  plot_data$Variable <- factor(plot_data$Variable, levels = affected)
  plot_data <- droplevels(plot_data)  # Drop unused levels
  plot_data <- plot_data[!is.na(plot_data$Estimate) & !is.na(plot_data$Variable), ]


  # Determine the range of event times
  min_event_time <- min(plot_data$EventTime)
  max_event_time <- max(plot_data$EventTime)

  # Create the plot
  ggplot(plot_data, aes(x = as.integer(EventTime), y = Estimate, color = Variable, group = Variable)) +
    geom_point(position = position_dodge(0.2), size = 1.5) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(0.2)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "solid", size = 0.5) +  # Add horizontal line at y=0
    scale_x_continuous(breaks = min_event_time:max_event_time) +  # Set x-axis breaks to integer values
    theme_minimal() +
    labs(x = "Event Time", y = "Estimate", color = "Variable") +
    theme(legend.position = "bottom")
}

# Example usage
# plotBLP(BLP_output, affected)



