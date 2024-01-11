## Declaring the MLDID class so that it creates
## an object to pass between the functions
#' MLDID Class
#'
#' This class is for storing the output of the main MLDID estimating procedure,
#' and is passed to follow-on functions. .
#'
#' @name MLDID
#' @import Matrix
#' @export
#' @slot attgt group-time ATTs
#' @slot cates Conditional Average Treatment Efffects
#' @slot scores double-robust scores
#' @slot gammas gammas used in the weighting procedure (Lu Nie Wager)
#' @slot positions positions of the variable (TO REMOVE)
#' @slot IDs Identifier of the observations
#' @slot params Parameters used by following functions
#' @slot processed_attgt Space for post-processed version of the estimated results
setClass("MLDID",
         slots = c(
           attgt = "list",
           cates = "dgCMatrix",  # assuming you use Matrix package for sparse matrices
           scores = "dgCMatrix",
           gammas = "dgCMatrix",
           positions = "dgCMatrix",
           IDs = "dgCMatrix",
           params = "list",
           processed_attgt = "list"  # New slot for processed data
         ))
