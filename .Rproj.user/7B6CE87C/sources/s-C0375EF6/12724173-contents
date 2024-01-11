### This file contains the Lu, Nie and Wager modified code 
## as function 
### To add: Regular DR ATT estimator (with scores extracted), and Zimmert 
## also added options for the ML models, and parameter tuning within the function#' 
#' ==========================================
#' @title Compute the Treatment Effects
#' 
#' @description This function runs machine learning for group-time ATTs (Average Treatment Effects on the Treated) 
#' and CATTs (Conditional Average Treatment Effects on the Treated) as implemented in the methodology 
#' of Lu, Nie, and Wager. It allows for various configurations and ML models to estimate treatment effects.
#'
#' @param X A data frame of covariates.
#' @param Y A numeric vector of outcomes.
#' @param Ti A numeric vector indicating time.
#' @param Si A numeric vector indicating subgroup membership.
#' @param constant_eff A character vector indicating whether the effect is constant or non-constant. 
#'        Accepted values are "constant" or "non_constant".
#' @param gamma Optional numeric value specifying the gamma parameter. If NULL, the default TMLE is used.
#' @param k_folds An integer indicating the number of folds for cross-validation. Default is 10.
#' @param lambda_choice A character string indicating the method for choosing lambda in the LASSO model. 
#'        Accepted values are "lambda.min" or "lambda.1se". Default is lambda.min
#' @param penalty_factor Optional numeric vector for penalty factors in the LASSO model.
#' @param tune_penalty A logical indicating whether to tune the penalty factors. Default is TRUE.
#' @param nu_model A character string specifying the model to use for nu estimation. 
#'        Accepted values are "rlearner" and "cf" (for causal forest).
#' @param sigma_model A character string specifying the model to use for sigma estimation. 
#'        Accepted values are "rlearner" and "cf" (for causal forest).
#' @param delta_model A character string specifying the model to use for delta estimation. 
#'        Accepted values are "glm" and "SL" (for SuperLearner).
#' @param t_func A logical indicating whether to use a functional form for t. Default is FALSE.
#'
#' @return A list containing the estimated treatment effects and associated statistics, 
#'         including TAU_hat, standard error estimates, tau_hat, score estimates, and y estimates.
#'
#' @examples
#' # Example usage:
#' # LNW_DiD(X = your_covariates, Y = your_outcomes, Ti = your_treatment_indicator, 
#' #         Si = your_subgroup_indicator)
#'
#' @export
LNW_DiD = function(X, Y, Ti, Si, constant_eff = c("constant","non_constant"), 
                   gamma = NULL, 
                   k_folds = 10,
                   lambda_choice = "lambda.min", 
                   penalty_factor = NULL, 
                   tune_penalty = TRUE, # default to tune the penalty factors for nu, sigma, delta   
                   nu_model = "rlearner",  # Default value, may also be "cf" for causal forest 
                   sigma_model = "rlearner",  # Default value, may also be "cf" for causal forest
                   delta_model = "glm", # Default, may also be "SL" for superlearner
                   t_func = FALSE) {
  # Check if nu_model is one of the allowed values
  if (!nu_model %in% c("cf", "rlearner")) {
    stop("Invalid value for nu_model. It must be either 'cf' or 'rlearner'.")
  }
  
  if (!sigma_model %in% c("cf", "rlearner")) {
    stop("Invalid value for nu_model. It must be either 'cf' or 'rlearner'.")
  }
  
  if (!delta_model %in% c("glm", "SL")) {
    stop("Invalid value for delta_model. It must be either 'glm' or 'CF'.")
  }
  
  constant_eff = match.arg(constant_eff)
  n = length(Ti)
  if (is.null(penalty_factor)){
    penalty_factor = rep(1, dim(X)[2])
  }
  
  # create foldid
  foldid = sample(rep(seq(k_folds), length = length(Ti)))
  
  # Fit values
  #print("Estimating nu")
  
  if (nu_model == "rlearner") {
    if (tune_penalty == FALSE) {
      nu_fit = rlearner::rlasso(X, Ti, Y, k_folds = k_folds, penalty_factor = penalty_factor)
      nu_hat = nu_fit$tau_hat
      t_hat = nu_fit$p_hat
      m_hat = nu_fit$m_hat
      print("nu is fitted with default rlearner.")
    } else {
      print("Tuning nu penalty factor ")
        penalty_factors <- c(0.01, 0.25, 0.5, 0.75, 0.99, 1)
        num_penalty_factors <- length(penalty_factors)
        cv_errors <- numeric(num_penalty_factors) 
        
        # Split data into k-folds
        folds <- cut(seq(1, nrow(X)), breaks=k_folds, labels=FALSE)
        
        # Perform cross-validation for each penalty factor
        for (i in seq_along(penalty_factors)) {
          fold_errors <- numeric(k_folds)
          
          for (k in 1:k_folds) {
            # Splitting the data: train and test
            test_indices <- which(folds == k)
            train_indices <- which(folds != k)
            
            X_train <- X[train_indices, ]
            Y_train <- Y[train_indices]
            Ti_train <- Ti[train_indices]
            
            X_test <- X[test_indices, ]
            Y_test <- Y[test_indices]
            Ti_test <- Ti[test_indices]
            
            # Fit the model on the training data
            fit <- rlearner::rlasso(X_train, Ti_train, Y_train, penalty_factor = rep(penalty_factors[i], ncol(X_train)))
            
            # Predict on the test data
            predictions <- predict(fit, newx = X_test)
            
            # Calculate error (e.g., Mean Squared Error)
            fold_errors[k] <- mean((Y_test - predictions)^2)
          }
          
          # Calculate the average error across all folds for the current penalty factor
          cv_errors[i] <- mean(fold_errors)
        }
        
        # Find the best penalty factor
        best_index <- which.min(cv_errors)
        best_penalty_factor <- penalty_factors[best_index]
      
        
        # Print the best alpha and penalty factor
        print(paste("Best nu penalty factor:", best_penalty_factor))
        print("Estimating nu with rlasso")
        nu_fit = rlearner::rlasso(X, Ti, Y, k_folds = k_folds, penalty_factor = rep(best_penalty_factor, ncol(X)))
        
        
        nu_hat = nu_fit$tau_hat
        t_hat = nu_fit$p_hat
        m_hat = nu_fit$m_hat
        print("nu is fitted.")
    }
    
  }
  else if (nu_model == "cf") {
    print("Estimating nu with causal forest")
    causal_forest_model = grf::causal_forest(X, Y, Ti, tune.parameters = c("mtry", "honesty.fraction", "imbalance.penalty"))
    nu_hat = causal_forest_model$predictions
    m_hat = causal_forest_model$Y.hat
    t_hat = causal_forest_model$W.hat
    
  }
  ## Replace t_hat with 0.5 since we split evenly by default (can change this as argument)
  #if (t_func == T) {
   # t_hat = nu_fit$p_hat
    #print(head(t_hat))
    #print(dim(t_hat))
  #}
  if (t_func == F) {
    print("Using constant T_hat")
    t_hat = rep(0.5, length(m_hat))
    }
   
  #print("Estimating sigma")
  if (sigma_model == "rlearner") {
    
    if (tune_penalty == FALSE) {
      penalty_factor = rep(1, dim(X)[2])
      sigma_fit = rlearner::rlasso(X, Si, Y, k_folds = k_folds, penalty_factor = penalty_factor)
      sigma_hat = sigma_fit$tau_hat
      s_hat = sigma_fit$p_hat
      print("sigma is fitted with default rlearner.")
      
    } else {
      print("Tuning sigma penalty factor ")
      penalty_factors <- c(0.01, 0.25, 0.5, 0.75, 0.99, 1)
      num_penalty_factors <- length(penalty_factors)
      cv_errors <- numeric(num_penalty_factors) 
      
      # Split data into k-folds
      folds <- cut(seq(1, nrow(X)), breaks=k_folds, labels=FALSE)
      
      # Perform cross-validation for each penalty factor
      for (i in seq_along(penalty_factors)) {
        fold_errors <- numeric(k_folds)
        
        for (k in 1:k_folds) {
          # Splitting the data: train and test
          test_indices <- which(folds == k)
          train_indices <- which(folds != k)
          
          X_train <- X[train_indices, ]
          Y_train <- Y[train_indices]
          Si_train <- Si[train_indices]
          
          X_test <- X[test_indices, ]
          Y_test <- Y[test_indices]
          Si_test <- Si[test_indices]
          
          # Fit the model on the training data
          fit <- rlearner::rlasso(X_train, Si_train, Y_train, penalty_factor = rep(penalty_factors[i], ncol(X_train)))
          
          # Predict on the test data
          predictions <- predict(fit, newx = X_test)
          
          # Calculate error (e.g., Mean Squared Error)
          fold_errors[k] <- mean((Y_test - predictions)^2)
        }
        
        # Calculate the average error across all folds for the current penalty factor
        cv_errors[i] <- mean(fold_errors)
      }
      
      # Find the best penalty factor
      best_index <- which.min(cv_errors)
      best_penalty_factor <- penalty_factors[best_index]
      
      
      # Print the best alpha and penalty factor
      print(paste("Best sigma penalty factor:", best_penalty_factor))
      print("Estimating sigma with rlasso")
      sigma_fit = rlearner::rlasso(X, Si, Y, k_folds = k_folds, penalty_factor = rep(best_penalty_factor, ncol(X)))
      
      
      sigma_hat = sigma_fit$tau_hat
      s_hat = sigma_fit$p_hat
      print("sigma is fitted.")
    }
    
  }
  else if (sigma_model == "cf"){
    print("Estimating sigma with causal forest")
    causal_forest_sigma = grf::causal_forest(X, Y, Si, tune.parameters = c("mtry", "honesty.fraction", "imbalance.penalty"))
    sigma_hat = causal_forest_sigma$predictions
    s_hat = causal_forest_sigma$W.hat
  }
  
  if (delta_model == "glm"){
    if (tune_penalty == FALSE) {
      penalty_factor = rep(1, dim(X)[2])
      delta_fit = glmnet::cv.glmnet(X, (Ti - t_hat) * (Si - s_hat), foldid = foldid, keep = TRUE, alpha = 1)
      delta_hat = as.matrix(delta_fit$fit.preval[,!is.na(colSums(delta_fit$fit.preval))])[, delta_fit$lambda == delta_fit$lambda.min]
      print("Delta fitted without tuning")
      } else {
        print("Tuning delta penalty factor ")
        
        # Define alpha values and penalty factors to test
        alpha_values <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
        penalty_factors <- c(0.01, 0.25, 0.5, 0.75, 0.99, 1)
        
        # Initialize variables to store results
        num_alphas <- length(alpha_values)
        num_penalty_factors <- length(penalty_factors)
        cv_errors <- matrix(NA, nrow = num_alphas, ncol = num_penalty_factors)
        models <- vector("list", length = num_alphas * num_penalty_factors)
        
        # Nested loop over alpha values and penalty factors
        for (i in seq_along(alpha_values)) {
          for (j in seq_along(penalty_factors)) {
            alpha <- alpha_values[i]
            penalty_factor <- rep(penalty_factors[j], ncol(X))
            model_index <- (i - 1) * num_penalty_factors + j
            
            # Fit the model
            fit <- glmnet::cv.glmnet(X, ((Ti - t_hat) * (Si - s_hat)), foldid = foldid, alpha = alpha, penalty.factor = penalty_factor, keep = TRUE)
            
            # Store the fitted model and its cross-validated error
            models[[model_index]] <- fit
            cv_errors[i, j] <- fit$cvm[fit$lambda == fit$lambda.min]
          }
        }
        
        # Identify the best combination of alpha and penalty factor
        best_indices <- which(cv_errors == min(cv_errors), arr.ind = TRUE)
        best_alpha <- alpha_values[best_indices[1, 1]]
        best_penalty_factor <- penalty_factors[best_indices[1, 2]]
        best_model_index <- (best_indices[1, 1] - 1) * num_penalty_factors + best_indices[1, 2]
        best_model <- models[[best_model_index]]
        
        # Print the best alpha and penalty factor
        print(paste("Best alpha:", best_alpha))
        print(paste("Best penalty factor:", best_penalty_factor))
        
        # Optionally refit the model using the best alpha and penalty factor on the entire dataset
        delta_fit <- glmnet::cv.glmnet(X, (Ti - t_hat) * (Si - s_hat), foldid = foldid, keep = TRUE, alpha = best_alpha)
        delta_hat = as.matrix(delta_fit$fit.preval[,!is.na(colSums(delta_fit$fit.preval))])[, delta_fit$lambda == delta_fit$lambda.min]
        
        print("delta is fitted with tuned glmnet.")
      }
  } else if (delta_model == "SL"){
    print("Estimating delta with SuperLearner")
    learners <- c("SL.glmnet", "SL.randomForest", "SL.xgboost", "SL.ranger", "SL.mean")    
    resp <- (Ti - t_hat) * (Si - s_hat)
    delta_fit = SuperLearner::CV.SuperLearner(Y = resp, X = data.frame(X), family = gaussian(),
                                cvControl = list(V = k_folds), innerCvControl = list(list(V=k_folds)),
                                SL.library = c("SL.glmnet", "SL.randomForest", "SL.xgboost")) #, "SL.ranger", "SL.mean"))
    delta_hat = delta_fit[["SL.predict"]]
    print("delta is fitted with SuperLearner.")
  }
  
  
  
  scaling = 1 - (delta_hat^2 / (s_hat *(1 - s_hat) * t_hat * (1 - t_hat)))
  A_hat = (Ti - t_hat - (delta_hat * (Si - s_hat)) / (s_hat * (Si - s_hat))) / scaling
  B_hat = (Si - s_hat - (delta_hat * (Ti - t_hat)) / (t_hat * (Ti - t_hat))) / scaling
  C_hat = Si * Ti - (s_hat + delta_hat / t_hat) * A_hat - (t_hat + delta_hat / s_hat) * B_hat - (s_hat * t_hat + delta_hat)
  
  S_hat = Y - m_hat - A_hat*nu_hat - B_hat*sigma_hat
  
  if (constant_eff == "constant") {
    input_data = as.data.frame(cbind(S_hat,C_hat))
    colnames(input_data) = c("S_hat", "C_hat")
    tau_fit_model = lm(S_hat~C_hat - 1, input_data)
    tau_fit_summary = summary(tau_fit_model)
    TAU_hat = tau_fit_summary$coefficients[1]
    std.err.est = sqrt(sandwich::vcovHC(tau_fit_model))
    
    return(list(TAU_hat = TAU_hat, std.err.est = std.err.est))
  } else if (constant_eff == "non_constant"){
    x_tilde = cbind(as.numeric(C_hat) * cbind(1, X))
    x_pred = cbind(1, X)
    n = length(Ti)
    k_folds = 10
    if (tune_penalty == FALSE) {
      p.fac = rep(1,dim(x_tilde)[2])
      p.fac[1] = 0
      foldid = sample(rep(seq(k_folds), length = length(Ti)))
      tau_fit = glmnet::cv.glmnet(x_tilde, S_hat, foldid = foldid, penalty.factor = p.fac)
      tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
      tau_hat = x_pred %*% tau_beta
    } else {
    x_pred = cbind(1, X)
    n = length(Ti)
    foldid = sample(rep(seq(k_folds), length = length(Ti)))
    ## going to try tuning the alpha 
    alphas2 = c(0, 0.01, 0.2, 0.4, 0.6, 0.8, 0.99, 1)
    penaltys2 <- c(0.01, 0.25, 0.5, 0.75, 0.99, 1)
    # Initialize variables to store results
    cv_errors2 <- matrix(NA, nrow = length(alphas2), ncol = length(penaltys2))
    models2 <- vector("list", length(alphas2) * length(penaltys2))
    print("Tuning tau parameters")
    # Nested loop over alpha values and penalty factors
    for (i in seq_along(alphas2)) {
      for (j in seq_along(penaltys2)) {
        alpha2 <- alphas2[i]
        penalty_factor2 <- rep(penaltys2[j], dim(x_tilde)[2])
        model_index2 <- (i - 1) * length(penaltys2) + j
        # Fit the model
        #print("Estimating tau")
        fit2 <- glmnet::cv.glmnet(x_tilde, S_hat, alpha = alpha2, foldid = foldid, penalty.factor = penalty_factor2)
        
        # Store the fitted model and its cross-validated error
        models2[[model_index2]] <- fit2
        cv_errors2[i, j] <- fit2$cvm[fit2$lambda == fit2$lambda.min]
      }
    }
    
    # Identify the best combination of alpha and penalty factor
    best_indices2 <- which(cv_errors2 == min(cv_errors2), arr.ind = TRUE)
    best_alpha2 <- alphas2[best_indices2[1, 1]]
    best_penalty_factor2 <- penaltys2[best_indices2[1, 2]]
    best_model_index2 <- (best_indices2[1, 1] - 1) * length(penaltys2) + best_indices2[1, 2]
    best_model2 <- models2[[best_model_index2]]
    
    # Print the best alpha and penalty factor
    print(paste("Best tau alpha:", best_alpha2))
    print(paste("Best tau penalty factor:", best_penalty_factor2))
    
    
    
    tau_fit = glmnet::cv.glmnet(x_tilde, S_hat, foldid = foldid, alpha = best_alpha2, penalty.factor = rep(best_penalty_factor2, dim(x_tilde)[2]))
    #tau_fit = CV.SuperLearner(Y = S_hat, X = data.frame(as.numeric(C_hat)*X), family = gaussian(),
     #                cvControl = list(V = k_folds), innerCvControl = list(list(V=k_folds)),
      #                SL.library = c("SL.glmnet", "SL.randomForest", "SL.xgboost"))
    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    #tau_beta = tau_fit[["SL.predict"]]
    tau_hat = x_pred %*% tau_beta
    #tau_hat = X %*% tau_beta
    print("tau is fitted with glmnet.")
    }
    if (is.null(gamma)){
      TAU_hat = mean(tau_hat)
      std.err.est = NA
    } else {
      y.est = m_hat + A_hat*nu_hat + B_hat*sigma_hat + C_hat*tau_hat
      TAU_hat = mean(tau_hat + gamma*(Y - y.est))
      score.est = tau_hat + gamma*(Y - y.est)
      std.err.est = sqrt((mean((tau_hat - TAU_hat)^2 + gamma^2 * (Y - y.est)^2))/n)
    }
    y.est = m_hat + A_hat*nu_hat + B_hat*sigma_hat + C_hat*tau_hat
    return(list(TAU_hat = TAU_hat, std.err.est = std.err.est,
                tau_hat = tau_hat, score.est = score.est, y.est = y.est))
  } else {
    stop("what effect?")
  }
  }
  
  
  