#' @title DBC criterion
#'
#' @description Deteminant-based criterion (DBC) is used to select the most
#'    appropriate working correlation struction in a generalized estimation
#'    equation (gee) setup.
#'
#' @param object \code{geeglm} object. The current implementation works only for
#'    "binomial" and "gaussian" families, and for "independence", "exchangeable",
#'     "ar1", "unstructured", and "toeplitz" correlation structures.
#'
#' @param corstr.toeplitz an indicator that takes  \code{TRUE} for
#'    "toeplitz" correlation structure
#'
#' @return It returns \code{dbc} criterion value
#'
#'  @importFrom geepack geeglm
#'
#' @seealso \code{geeglm}
#'
#' @references Jaman A, Latif AHMM, Bari W, and Wahed A (2016). A determinant‐based
#'    criterion for working correlation structure selection in generalized
#'    estimating equations. Statistics in Medicine 35(11), 1819-1833.
#'
#' @export dbc
#'
#' @examples
#'
#' library(geepack)
#' data(ohio)
#' fit <- geeglm(resp ~ age + smoke, id=id, data=ohio,
#'         family=binomial, corstr="exch")
#' dbc(fit)
#'
#'


#############################################################################
## DESCRIPTION:
##
## One needs to source all the following functions for calculating DBC
## for a given GEE model. To obtain the value of DBC criterion we need
## to run DBC.fun function with a geeglm object as argument. The current
## function only works for "binomial" and "gaussian" families and for
## "independence", "exchangeable", "ar1", "unstructured", and "toeplitz"
## correlation structures. If someone wants to obtain DBC for toeplitz
## correlation structure, then he/she must define corstr.toeplitz = TRUE
## in DBC.fun and should provide the geeglm object, which has been fitted
## with toeplitz structure.
##
## EXAMPLE:
##
## myModelFit <- geeglm(formula = myFormula, data = myData, id = myData$id,
##                  family = ?binomial?, corstr = ?ar1?))
## DBC <- DBC.fun(object = myModelFit)
##
#############################################################################

dbc <- function(object, corstr.toeplitz = FALSE) {
  #
   if (summary(object)$error == 1) {
     message("geeglm fit did not attain convergence (summary.geeglm$error = 1)")
     return(NA)
   }
   family <- object$family$family
   if (!(family %in% c("gaussian", "binomial"))) {
      message("Family in geeglm must be gaussian or binomial")
     return(NA)
   }
   data.st <- data.frame(id = object$id,
     stats::model.frame(formula = object$formula,
         data = object$data))
   N <- length(unique(data.st$id))
   Ni <- table(data.st$id)
   fit <- summary(object)
   b <- fit$coef[,1]
   phi <- unlist(fit$dispersion[1])
   alpha <- unlist(fit$corr[1])
   if (corstr.toeplitz == TRUE && length(alpha) != max(Ni) - 1) {
      message("When corstr.toeplitz is TRUE geeglm object must have toeplitz structure")
     return(NA)
   }
   if (corstr.toeplitz == TRUE) {
      corsr = "toeplitz"
   } else {
      corsr = object$corstr
   }
   Vm <- object$geese$vbeta.naiv
   Vm_inv <- solve(Vm)
   Vs_bc <- Vs.biasC(V.model = Vm, data = data.st, N = N, Ni = Ni, beta = b,
         alpha = alpha, phi = phi, corsr = corsr, family = family)
   dbc <- det(Vm_inv %*% Vs_bc)
   names(dbc) <- "dbc"
   return(dbc)
}
