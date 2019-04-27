
#############################################################################
##
## Vs.biasC function calculates the bias corrected sandwich estimator.
##
#############################################################################

Vs.biasC <- function(V.model, data, N, Ni, beta, alpha, phi, corsr, family){
   dat <- split(data, data[,1])
   sum_det <- 0
   if (corsr == "unstructured" || corsr == "toeplitz") {
      wcm.hat = corMatF(alpha = alpha, ni = max(Ni), corstr = corsr)
   }
   for (i in 1:N) {
      y <- as.matrix(dat[[i]][,2])
      x <- as.matrix(cbind(1,dat[[i]][,3:length(dat[[1]])]))
      mn <- mean.fun(x = x, beta = beta, family = family)
      if (corsr == "unstructured" || corsr == "toeplitz") {
         wcm.hat.i = wcm.hat[1:Ni[i],1:Ni[i]]
      }
      alpha.i <- switch(corsr,
         "independence" = NULL,
         "exchangeable" = alpha,
         "ar1" = alpha,
         "unstructured" = wcm.hat.i[lower.tri(wcm.hat.i)]
      )
      if (corsr == "toeplitz") {
         if (Ni[i] == 1) {
            alpha.i = NULL
         } else {
            alpha.i = wcm.hat.i[1,2:Ni[i]]
         }
      }
      ei <- y - mn
      ei_tr <- t(ei)
      Di <- calc.D(nrc = Ni[i], x = x, beta = beta, family = family)
      Di_tr <- t(Di)
      Vi <- calc.V(alfa = alpha.i, fi = phi, corsr = corsr, nrc = Ni[i],
      mean.vect = mn, family = family)
      Vi_inv <- solve(Vi)
      Hi <- Di %*% V.model %*% Di_tr %*% solve(Vi)
      I <- diag(Ni[i])
      I_Hi_inv <- solve(I - Hi)
      sum_det <- sum_det + Di_tr %*% Vi_inv %*% I_Hi_inv %*% ei %*% ei_tr %*%
                   I_Hi_inv %*% Vi_inv %*% Di
   }
   Vs.biasC <- V.model %*% sum_det %*% V.model
   return(Vs.biasC)
}

corMatF <- function(alpha, ni, corstr){
   switch(corstr,
      "exchangeable" = stats::toeplitz(c(1, rep(alpha, ni - 1))),
      "ar1" = stats::toeplitz(alpha^(0:(ni - 1))),
      "unstructured" = vec2uMat(alpha, ni),
      "independence" = diag(ni),
      "toeplitz" = stats::toeplitz(c(1, alpha))
   )
}

vec2uMat <- function(alpha, ni){
   x <- matrix(1, ni, ni)
   x[lower.tri(x)] <- alpha
   x[upper.tri(x)] <- t(x)[upper.tri(x)]
   return(x)
}

mean.fun <- function(x, beta, family){
   switch(family,
      "binomial" = clogit(x %*% beta),
      "gaussian" = x %*% beta
   )
}

clogit <- function(x) exp(x)/(1 + exp(x))

calc.D <- function(nrc, x, beta, family){
   if (family == "gaussian") {
      D = x
   }
   if (family == "binomial") {
      l.pred <- x %*% beta
      c <- exp(l.pred)/(1 + exp(l.pred))^2
      d <- matrix(c, ncol = length(beta), nrow = nrc)
      D <- d*x
   }
   return(D)
}

calc.V <- function(alfa, fi, corsr, nrc, mean.vect, family){
   R <- corMatF(alpha = alfa, ni = nrc, corstr = corsr)
   Ai <- diag(c(fi*var.fun(mean = mean.vect, family = family)), nrc)
   V <- mat.power(Ai, 1/2) %*% R %*% mat.power(Ai, 1/2)
   return(V)
}


var.fun <- function(mean, family){
   switch(family,
     "gaussian" = rep(1, length(mean)),
     "binomial" = mean*(1 - mean)
   )
}

mat.power <- function(mat, power){
   e.vect <- eigen(mat)$vectors
   e.val.p <- (eigen(mat)$values)^power
   res.mat <- e.vect %*% diag(e.val.p, length(e.val.p)) %*% t(e.vect)
   return(res.mat)
}

