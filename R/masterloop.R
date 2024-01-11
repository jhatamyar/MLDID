#'
#' ==========================================
#' @title Compute Treatment Effect Parameters
#'
#' @description runs the machine learning for the group-time ATTs and CATTs. See README for example usage.
#'
#' @import methods
#' @import BMisc
#'
#' @param outcome Name of the outcome variable in the data.
#' @param group Name of the group variable in the data, typically indicating treatment assignment.
#' @param time Name of the time variable in the data.
#' @param id_name Name of the identifier variable in the data.
#' @param data The dataframe to be used for the analysis.
#' @param xformla A formula object specifying the covariates to be used in the model.
#' @param tune_penalty Logical flag indicating whether to tune the penalty parameters (TRUE) or not (FALSE).
#' Parameters to be tuned are in nu, sigma, and delta - see the function LNW_DID for details
#' @param nu_model Model to be used for the nu component. Possible values are "rlearner" (default) or "cf" for causal forest.
#' @param sigma_model Model to be used for the sigma component. Possible values are "rlearner" (default) or "cf" for causal forest.
#' @param delta_model Model to be used for the delta component. Possible values are "glm" (default) or "SL" for superlearner.
#' @param t_func A function specifying whether to estimate t_hat. Default is (TRUE), if (FALSE), will replace with 0.5, since there is always one observation pre/post
#' NOTE that you MUST specify this as FALSE if there are no time-varying covariates in the data
#' @param verbose Logical flag indicating whether to print detailed information during function execution (default is TRUE).
#' @return An object of class `MLDID` containing the estimated treatment effects and other
#'         relevant information. This object includes the following components:
#'         - attgt: A list of estimated treatment effects for each group and time period, using both the MLDID
#'         and DRDID (Callaway & Sant'Anna) estimators. The MLDID estimates are 'TAU-hat' and the DRDID are 'attgt'.
#'         - cates: A matrix of Conditional Average Treatment Effects on the Treated.
#'         - scores: A matrix of double-robust scores for the Treated.
#'         - gammas: A matrix of estimated gammas
#'         - positions: A matrix indicating the positions of each unit in the data.
#'         - IDs: A matrix of IDs for each unit in the data.
#'         - params: A list of parameters used in the simulation, including data, time lists,
#'                  group lists, and number of periods, groups, and units.
#' @export
MLDID <- function(outcome,
                  group,
                  time,
                  id_name,
                  data,
                  xformla = xformla,
                  tune_penalty = TRUE,
                  nu_model = "rlearner", # Default, can also be "cf" for causal forest
                  sigma_model = "rlearner", #Default, can also be "cf" for causal forest
                  delta_model = "glm", # Default, can also be "SL" for superlearner
                  t_func = TRUE,
                  verbose = TRUE) {
  # call the pre-processing function from _helperfuncs.R
  ret=did_datapreprocess(outcome = outcome,
                        group = group,
                        time = time,
                        id_name = id_name,
                        data = data)

  ## unpack params
  dta   <- ret$data
  tlist <- ret$tlist
  glist <- ret$glist
  nT    <- ret$nT
  nG    <- ret$nG
  n     <- ret$n

  # create place holder in lists
  counter <- 1
  attgt.list <- list()

  # number of time periods
  tlist.length <- length(tlist)
  tfac <- 0 ## corresponds to "universal" base period

  ## also need list of event times for later use
  temp.data <- dta[dta$G != 0,]
  eseq <- unique(temp.data$period - temp.data$G)
  e.times <- eseq[order(eseq)]

  # 3-dimensional array which will store cates
  # across groups and times

  ## Note: can add to this to store score fn once I figure out what it is
  cates <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)
  scores <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)
  y.hats <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)
  gammas <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)

  ## for keeping track later
  IDs <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)
  positions <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)

  # loop over groups
  for (g in 1:(nG)) { ## must subtract one here because there is no year at the end without TX?

    ## set up G once: will be a 1 if this is the current group
    dta$disG <- 1*(dta$G == glist[g])
    # loop over time periods
    for (t in 1:tlist.length) {


      # use same base period as for post-treatment periods
      pret <- tail(which( (tlist) < glist[g]),1)


      # Print statements wrapped in if (verbose)
      if (verbose) {
        cat(paste("current period:", tlist[(t)]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pretreatment period to be", tlist[pret]), "\n")
      }


      # use "not yet treated as control"
      dta$disC <- 1 * ((dta$G == 0) |
                         ((dta$G > (tlist[max(t,pret)])) &
                            (dta$G != glist[g])))
      # check if in post-treatment period
      if ((glist[g]<=tlist[(t)])) {

        # update pre-period if in post-treatment period to
        # be  period (g-delta-1)
        pret <- tail(which( (tlist) < glist[g]),1)

        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, code will
          # jump out of this loop
          break
        }
      }

      # if we are in period (g-1), normalize results to be equal to 0
      # and break without computing anything
      if (tlist[pret] == tlist[(t)]) {
        ### NEED TO UPDATE THIS
        attgt.list[[counter]] <- list(TAU_hat = 0,
                                      attgt = 0,
                                      taus = 0,
                                      std.err.est = 0,
                                      group=glist[g],
                                      year=tlist[(t)],
                                      post=0)
        #inffunc[,counter] <- rep(0,n)
        counter <- counter+1
        next
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------
      ###--------------------------------
      # note we need to create Ti, Si, Y
      ###--------------------------------

      # post treatment dummy variable
      post.treat <- 1*(glist[g]<=tlist[t])

      # get dataset with current period and pre-treatment period
      # total number of units (not just included in G or C)
      #disdat <- data[data[,tname] == tlist[t+tfac] | data[,tname] == tlist[pret],]
      disdat <- dta[(dta$period==tlist[t] | dta$period==tlist[pret]),]
      ## disdat2 is in the form for DRDID
      disdat.SA <- panel2cs2(disdat, "Y", "id", "period", balance_panel=FALSE)

      # still total number of units (not just included in G or C)
      n <- nrow(disdat)
      n.SA <- nrow(disdat.SA)

      # pick up the indices for units that will be used to compute ATT(g,t)
      disidx <- disdat$disG==1 | disdat$disC==1
      disidx.SA <- disdat.SA$disG==1 | disdat.SA$disC==1

      # pick up the data that will be used to compute ATT(g,t)
      disdat <- disdat[disidx,]
      disdat.SA <- disdat.SA[disidx.SA,]

      n1 <- nrow(disdat) # num obs. for computing ATT(g,t)
      n1.SA <- nrow(disdat.SA)

      # drop missing factors MAYBE DELETE
      disdat <- droplevels(disdat)
      disdat.SA <- droplevels(disdat.SA)

      # give short names for data in this iteration
      G <- disdat$disG
      C <- disdat$disC
      G.SA <- disdat.SA$disG
      C.SA <- disdat.SA$disC

      ## trying the new code version
      #
      Ypre <- if(tlist[(t)] > tlist[pret]) disdat.SA$.y0 else disdat.SA$.y1
      Ypost <- if(tlist[(t)] > tlist[pret]) disdat.SA$.y1 else disdat.SA$.y0


      # matrix of covariates
      covariates <- model.matrix(xformla, data=disdat)
      # remove column of 1s due to difference in package
      covariates <- covariates[,2:ncol(covariates)]

      covariates.SA <- model.matrix(xformla, data=disdat.SA)
      ##

      Y <- disdat$Y
      #Si: INDICATES IF IN TREATED GROUP
      Si <- disdat$disG
      ## timing indicator should be indicating if in second period
      disdat$disT <- as.numeric(disdat$period == tlist[t])
      Ti <- disdat$disT


      #-----------------------------------------------------------------------------
      # code for actually computing att(g,t)
      #-----------------------------------------------------------------------------

      ### TO DO: give option to change gamma
      ### trying to estimate a gamma for weights
      ### COPYING the minimax solver by hand, for some reason the code doesn't work
      #nobs = nrow(X)
      nobs = nrow(covariates)
      #pobs = ncol(X)
      pobs = ncol(covariates)
      gg = CVXR::Variable(nobs + 4)

      zeta=0.5
      objective = (1 - zeta) * sum(gg[1:nobs]^2) + zeta * sum(gg[nobs + 1:4]^2)


      contraints = list(
        sum(gg[1:nobs]) == 0,
        t(covariates) %*% gg[1:nobs] <= gg[nobs + 1],
        -t(covariates) %*% gg[1:nobs] <= gg[nobs + 1],
        t(covariates) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
        -t(covariates) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
        t(covariates) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
        -t(covariates) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
        sum(Ti * Si * gg[1:nobs]) == 1,
        t(covariates) %*% (Ti * Si * gg[1:nobs]) <= colMeans(covariates) + gg[nobs + 4],
        - t(covariates) %*% (Ti * Si * gg[1:nobs]) <= - colMeans(covariates) + gg[nobs + 4]
      )

      cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
      cvx.output = solve(cvx.problem, solver = "ECOS", verbose = F)
      result = cvx.output$getValue(gg)
      gamma = nobs * result[1:nobs]

      TR_fit <- LNW_DiD(X = covariates,
                                Y = Y,
                                Ti = Ti,
                                Si = Si,
                                gamma = gamma,
                                lambda_choice = "lambda.1se",
                                constant_eff = "non_constant",
                                k_folds = 5,
                                t_func = t_func,
                                tune_penalty = tune_penalty,
                                nu_model = nu_model,
                                sigma_model = sigma_model,
                                delta_model = delta_model
                        )

      attgt <- DRDID::drdid_panel(Ypost, Ypre, G.SA,
                                  covariates = covariates.SA,
                                  #i.weights=w,
                                  boot=FALSE, inffunc=TRUE
      )

      # save results for this att(g,t)
      attgt.list[[counter]] <- list(TAU_hat = TR_fit$TAU_hat,
                                    attgt = attgt$ATT,
                                    attgt.se = attgt$se,
                                    taus = TR_fit$tau_hat,
                                    std.err.est = TR_fit$std.err.est,
                                    group=glist[g],
                                    year=tlist[(t)],
                                    post=post.treat)

      ## save cates
      cate.temp <- rep(NA, n1)
      cate.temp[disidx] <- TR_fit$tau_hat
      cates[,counter] <- cate.temp

      ## save score.est
      score.temp <- rep(NA, n1)
      score.temp[disidx] <- TR_fit$score.est
      scores[,counter] <- score.temp

      ## save y.est
      y.temp <- rep(NA, n1)
      y.temp[disidx] <- TR_fit$y.est
      y.hats[,counter] <- y.temp

      # save gamma
      g.temp <- rep(NA, n1)
      g.temp[disidx] <- gamma
      gammas[,counter] <- g.temp

      ## save positions and IDs for later
      positions.temp <- rep(NA, n1)
      positions.temp[disidx] <- disdat$period
      positions[,counter] <- positions.temp

      IDs.temp <- rep(NA, n1)
      IDs.temp[disidx] <- disdat$id
      IDs[,counter] <- IDs.temp

      # update counter
      counter <- counter+1
    } # end looping over
  } # end looping over g

  return(new("MLDID",
             attgt = attgt.list,
             cates = cates,
             scores = scores,
             gammas = gammas,
             positions = positions,
             IDs = IDs,
             params = list("data" = dta,
                           "tlist" = tlist,
                           "glist" = glist,
                           "nT" = nT,
                           "nG" = nG,
                           "n" = n,
                           "e.times" = e.times)))

}






#'@title Obtaining oracle effects for simulations
#'
#'@description This obtains the true ATTs and CATTs by direct calculation, using the known values
#'of the potential outcomes under treatment and control.
#'
#'
#' @param dta A data frame containing the simulation data. This should include variables for
#'            treatment status, subgroup membership, and potential outcomes under treatment and control.
#' @param xformla A formula or other specification of the model used in the simulation,
#'                indicating how the covariates are related to the outcomes.
#'
#' @return An object of class "MLDID" containing the estimated treatment effects and other
#'         relevant information. This object includes the following components:
#'         - attgt: A list of estimated treatment effects for each group and time period.
#'         - cates: A matrix of Conditional Average Treatment Effects on the Treated.
#'         - positions: A matrix indicating the positions of each unit in the data.
#'         - IDs: A matrix of IDs for each unit in the data.
#'         - params: A list of parameters used in the simulation, including data, time lists,
#'                  group lists, and number of periods, groups, and units.
#'
#' @examples
#' # Example usage:
#' # didMLloop_oracle(dta = your_data, xformla = your_formula)
#'
#' @import methods
#'
#' @export
didMLloop_oracle <- function(dta = dta,
                      xformla = xformla) {

  # call the pre-processing function from _helperfuncs.R
  ret=did_datapreprocess_sim(dta)

  ## unpack params
  dta   <- ret$data
  tlist <- ret$tlist
  glist <- ret$glist
  nT    <- ret$nT
  nG    <- ret$nG
  n     <- ret$n

  # place holder in lists
  counter <- 1
  attgt.or.list <- list()

  # number of time periods
  tlist.length <- length(tlist)
  tfac <- 0 ## corresponds to "universal" base period

  ## also need list of event times for later use
  temp.data <- dta[dta$G != 0,]
  eseq <- unique(temp.data$period - temp.data$G)
  e.times <- eseq[order(eseq)]

  # 3-dimensional array which will store cates
  # across groups and times

  ## Note: can add to this to store score fn once I figure out what it is
  cates.or <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)

  ## for keeping track later
  IDs <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)
  positions <- Matrix::Matrix(data=0,nrow=n*2, ncol=(nG)*(nT), sparse=TRUE)

  for (g in 1:(nG)) { ## must subtract one here because there is no year at the end without TX?

    ## set up G once: will be a 1 if this is the current group
    dta$disG <- 1*(dta$G == glist[g])
    # loop over time periods
    for (t in 1:tlist.length) {


      # use same base period as for post-treatment periods
      pret <- tail(which( (tlist) < glist[g]),1)


      # print the details of which iteration we are on

      cat(paste("current period:", tlist[(t)]), "\n")
      cat(paste("current group:", glist[g]), "\n")
      cat(paste("set pretreatment period to be", tlist[pret]), "\n")

      # use "not yet treated as control"
      dta$disC <- 1 * ((dta$G == 0) |
                          ((dta$G > (tlist[max(t,pret)])) &
                             (dta$G != glist[g])))
      # check if in post-treatment period
      if ((glist[g]<=tlist[(t)])) {

        # update pre-period if in post-treatment period to
        # be  period (g-delta-1)
        pret <- tail(which( (tlist) < glist[g]),1)

        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, code will
          # jump out of this loop
          break
        }
      }

      # if we are in period (g-1), normalize results to be equal to 0
      # and break without computing anything
      if (tlist[pret] == tlist[(t)]) {

        attgt.or.list[[counter]] <- list(TAU.or = 0,
                                         taus.or = 0,
                                         group=glist[g],
                                         year=tlist[(t)],
                                         post=0)

        counter <- counter+1
        next
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------
      ###--------------------------------
      # note we need to create Ti, Si, Y
      ###--------------------------------

      # post treatment dummy variable
      post.treat <- 1*(glist[g]<=tlist[t])

      # get dataset with current period and pre-treatment period
      # total number of units (not just included in G or C)
      disdat <- dta[(dta$period==tlist[t] | dta$period==tlist[pret]),]

      # still total number of units (not just included in G or C)
      n <- nrow(disdat)

      # pick up the indices for units that will be used to compute ATT(g,t)
      disidx <- disdat$disG==1 | disdat$disC==1

      # pick up the data that will be used to compute ATT(g,t)
      disdat <- disdat[disidx,]

      n1 <- nrow(disdat) # num obs. for computing ATT(g,t)

      # drop missing factors MAYBE DELETE
      disdat <- droplevels(disdat)

      # give short names for data in this iteration
      G <- disdat$disG
      C <- disdat$disC

      Y <- disdat$Y
      y1 <- disdat$y1
      y0 <- disdat$y0

      Si <- disdat$disG
      ## timing indicator should be indicating if in second period
      disdat$disT <- as.numeric(disdat$period == tlist[t])
      Ti <- disdat$disT


      #-----------------------------------------------------------------------------
      # code for actually computing att(g,t)
      #-----------------------------------------------------------------------------

      att.or <- mean(y1[Ti==1 & G == 1] - y0[Ti==1 & G == 1])

      disdat$cate.or <- ifelse((disdat$disT == 1 & disdat$disG == 1), (disdat$y1 - disdat$y0), NA)
      #disdat$cate.or <- disdat$y1 - disdat$y0
        #y1[Ti==1 & G == 1] - y0[Ti==1 & G == 1]
      # pick up the indices for units that will be used to compute cates
      #discateidx <- disdat$disT==1 & disdat$disG==1
      # pick up the data that will be used to compute cates
      #discatedat <- disdat[discateidx,]

      #ncates <- nrow(discatedat) # num obs. for computing ATT(g,t)

      # pick up the data that will be used to compute ATT(g,t)
      #disdat <- disdat[disidx,]
      #cate.or <- y1 - y0

      # save results for this att(g,t)
      attgt.or.list[[counter]] <- list(TAU.or = att.or,
                                       taus.or = disdat$cate.or,
                                       group=glist[g],
                                       year=tlist[(t)],
                                       post=post.treat)

      ## save cates
      cate.temp <- rep(NA, n1)
      cate.temp[disidx] <- disdat$cate.or
      cates.or[,counter] <- cate.temp


      ## save positions and IDs for later
      positions.temp <- rep(NA, n1)
      positions.temp[disidx] <- disdat$period
      positions[,counter] <- positions.temp

      IDs.temp <- rep(NA, n1)
      IDs.temp[disidx] <- disdat$id
      IDs[,counter] <- IDs.temp

      # update counter
      counter <- counter+1
    } # end looping over
  } # end looping over g

  return(new("MLDID",
             attgt = attgt.or.list,
             cates = cates.or,
             positions = positions,
             IDs = IDs,
             params = list("data" = dta,
                           "tlist" = tlist,
                           "glist" = glist,
                           "nT" = nT,
                           "nG" = nG,
                           "n" = n,
                           "e.times" = e.times)))


}
