#' @title DBC criterion
#'
#' @description Calculates deteminant-based criterion (DBC), which can be 
#'    used to select the most appropriate working correlation struction 
#'    in generalized estimation equation (gee). Either the arguments of 
#'    \code{geeglm} or its output can be used as arguments of \code{dbc}.
#'
#' @param object \code{geeglm} output object and the current implementation 
#'     works only for "binomial" and "gaussian" families, and for 
#'     "independence", "exchangeable", "ar1", "unstructured", and 
#'     "toeplitz" correlation structures.
#' 
#' @param formula see corresponding documentation of \code{glm}
#' 
#' @param family see corresponding documentation of \code{glm}
#' 
#' @param data see corresponding documentation of \code{glm}
#' 
#' @param id see corresponding documentation of \code{geeglm}
#' 
#' @param corstr see corresponding documentation of \code{geeglm}
#' 
#' @param corstr.toeplitz an indicator that takes  \code{TRUE} for
#'    "toeplitz" correlation structure
#'
#' @param ... further arguments of \code{geeglm} to be passed 
#' 
#' @return rerurns \code{dbc} criterion value and the output of \code{geeglm}
#'
#' @importFrom geepack geeglm
#'
#' @seealso \code{geeglm}, \code{glm}
#'
#' @references Jaman A, Latif AHMM, Bari W, and Wahed A (2016). A determinant‐based
#'    criterion for working correlation structure selection in generalized
#'    estimating equations. Statistics in Medicine 35(11), 1819-1833.
#'
#' @export dbc
#'
#' @examples
#' library(geepack)
#' data(ohio)
#' #
#' fit <- geeglm(resp ~ age + smoke, id=id, data=ohio,
#'         family=binomial, corstr="exch")
#' dbc(fit)
#' #
#' # geeglm arguments can also be used as arguments 
#'
#' fit2 <- dbc(formula = resp ~ age + smoke, id=id, data=ohio,
#'         family=binomial, corstr="exch")
#'         
 dbc <- function(object = NULL, formula, family, data, id, 
   corstr, corstr.toeplitz = FALSE, ...) {
   #
   if (is.null(object)) {
     #formula <- as.formula(formula)
     object <- geepack::geeglm(formula = formula, family = family, id = id, 
       corstr = corstr,
       data = data, ...)
   }
   if (summary(object)$error == 1) {
     message(
       "geeglm fit did not attain convergence (summary.geeglm$error = 1)"
       )
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
      message(
        "When corstr.toeplitz is TRUE geeglm object must have toeplitz structure"
        )
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
   #names(dbc) <- "dbc"
   return(list(fit = object, dbc = dbc))
}
