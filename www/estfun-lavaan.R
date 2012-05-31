estfun.lavaan <- function(x)
{
  ## observed data
  X <- x@Data@X
  ## number variables/sample size
  nvar <- x@Model@nvar
  nobs <- x@SampleStats@ntotal
  ntab <- unlist(x@SampleStats@nobs)

  ## number of groups
  ngrp <- x@Model@ngroups

  ## Define matrix that will hold all scores
  ## Assumes means are estimated ## FIXME: check and otherwise throw error
  scores.H0 <- matrix(NA, nobs, x@Model@nx.free)
  
  for (i in 1:ngrp) {
    ## fitted moments
    if(ngrp > 1){
      moments <- fitted(x)[[i]]
    }else{
      moments <- fitted(x)
    }
    Sigma.hat <- moments$cov
    ## deal with wishart likelihood? FIXME: check
    if (x@Options$likelihood == "wishart"){
      Sigma.hat <- ((ntab[i] - 1)/(ntab[i])) * Sigma.hat
    }
    Mu.hat <- moments$mean
    Sigma.inv <- solve(Sigma.hat)

    ## Junk matrices for multiplication
    J <- matrix(1, 1L, ntab[i]) ## FIXME: needed? better maybe rowSums/colSums?
    J2 <- matrix(1, nvar[i], nvar[i])
    diag(J2) <- 0.5

    ## scores.H1 (H1 = saturated model)
    mean.diff <- t(t(X[[i]]) - Mu.hat %*% J)

    dx.Mu <- -1 * mean.diff %*% Sigma.inv

    dx.Sigma <- t(apply(mean.diff, 1L,
      function(x) lavaan:::vech(- J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat) %*% Sigma.inv))))

    scores.H1 <- cbind(dx.Mu, dx.Sigma)
    
    ## scores.H0
    Delta <- lavaan:::computeDelta(x@Model)[[i]]

    ## assign to rows of overall score matrix
    wi <- x@Data@case.idx[[i]]
    scores.H0[wi, ] <- scores.H1 %*% Delta
  }

  ## return gradient of log-likelihood
  rval <- -scores.H0
  colnames(rval) <- names(coef(x))
  return(rval)
}
