################################################################
## A, B as defined in Vuong Eq (2.1) and (2.2)
################################################################
calcAB <- function(object){
  ## Eq (2.1)
  samplestats <- object@SampleStats
  ntot <- samplestats@ntotal

  A <- solve(ntot * vcov(object))

  ## Eq (2.2)
  tmp <- estfun.lavaan(object)
  tmp2 <- crossprod(tmp)/ntot
  B <- matrix(tmp2, nrow(A), nrow(A))

  list(A=A, B=B)
}

## And a function to get the cross-product from Eq (2.7)
calcBcross <- function(object1, object2){
  ## Get Eq (2.7)
  tmp1 <- estfun.lavaan(object1)
  tmp2 <- estfun.lavaan(object2)

  crossprod(tmp1, tmp2)/object1@SampleStats@ntotal
}


################################################################
## Calculating W, Vuong Eq (3.6)
################################################################
calcLambda <- function(object1, object2) {
    AB1 <- calcAB(object1)
    AB2 <- calcAB(object2)
    Bc <- calcBcross(object1, object2)

    W <- cbind(rbind(-AB1$B %*% solve(AB1$A), t(Bc) %*% solve(AB1$A)),
               rbind(-Bc %*% solve(AB2$A), AB2$B %*% solve(AB2$A)))

    lamstar <- eigen(W, only.values=TRUE)$values
    ## Discard imaginary part, as it only occurs for tiny eigenvalues?
    # as.numeric(lamstar)^2
    ## using Mod
    Mod(lamstar)^2
}


################################################################
## Getting log-likelihood for individual cases
################################################################

get.ll <- function(object){
  samplestats <- object@SampleStats
  llvec <- rep(NA, samplestats@ntotal)

  for(g in 1:samplestats@ngroups) {
    if (samplestats@ngroups > 1){
      moments <- fitted(object)[[g]]
    } else {
      moments <- fitted(object)
    }
    Sigma.hat <- unclass(moments$cov)
    ## To ensure it is symmetric; needed?
    Sigma.hat <- (Sigma.hat + t(Sigma.hat))/2
    ## Which cases correspond to this group?
    grpind <- object@Data@case.idx[[g]]
    
    if(!samplestats@missing.flag) { # complete data
      if(object@Model@meanstructure) { # mean structure
        Mu.hat <- unclass(moments$mean)
      } else {
        ## set mean structure to sample estimates
        Mu.hat <- apply(object@Data@X[[g]], 2, mean)
      }
      llvec[grpind] <- dmnorm(object@Data@X[[g]], Mu.hat, Sigma.hat, log=TRUE)
    } else { # incomplete data
      nsub <- ntab[g]
      M <- samplestats@missing[[g]]
      MP1 <- object@Data@Mp[[g]]
      pat.idx <- match(MP1$id, MP1$order)
      tmpll <- rep(NA, nsub)

      Mu.hat <- unclass(moments$mean)
      nvar <- ncol(samplestats@cov[[g]])

      for(p in 1:length(M)) {
        ## Data
        X <- M[[p]][["X"]]
        var.idx <- M[[p]][["var.idx"]]

        tmpll[pat.idx==p] <- dmnorm(X, Mu.hat[var.idx], Sigma.hat[var.idx, var.idx], log=TRUE)
      }

      llvec[grpind] <- tmpll
    } # incomplete
  } # group
  llvec
}

################################################################
## Calculate and test \hat{\omega}^2; Eq (4.2)
################################################################

vuong <- function(object1, object2, data) {
    llA <- get.ll(object1)
    llB <- get.ll(object2)

    ## Eq (4.2)
    n <- nrow(data)
    omega.hat.2 <- (n-1)/n * var(llA - llB)

    ## Get p-value of weighted chi-square dist
    ## Need to install the dr package
    lamstar2 <- calcLambda(object1, object2)
    tmp <- dr.pvalue(lamstar2, n * omega.hat.2)
    pOmega <- tmp[[4]]

    ## Calculate and test LRT; Eq (6.4)
    lr <- sum(llA - llB)
    teststat <- (1/sqrt(n)) * lr/sqrt(omega.hat.2)

    ## Not needed: We know it will prefer the better fitting model:
    ## sign of teststat: 1 for object1, -1 for object2
    pref <- sign(teststat)

    ## 1-tailed p-value, because the better fitting model can only
    ## be preferred
    pLRT <- pnorm(abs(teststat), lower.tail=FALSE)

    ## Interval for BIC differences
    bicA <- BIC(object1)
    bicB <- BIC(object2)
    bicdiff <- bicA - bicB

    ## 90% two-sided CI to match one-sided test with alpha=.05
    ci <- bicdiff + qnorm(c(.05,.95))*sqrt(n * 4 * omega.hat.2)
    ## I believe we don't need this, because the null hypothesis is
    ## that BIC difference=0, not that LR=0:
    ## Adjust away the difference in parameters, so can compare to 0?
    ##k <- length(coef(object1))
    ##q <- length(coef(object2))
    ##adjci <- ci - (k - q)*log(n)

    list(omega = omega.hat.2, p_Omega=pOmega, z_LRT = teststat, p_LRT=pLRT, pref=pref, BICint=ci) #, BICintadj=adjci)
}

## check that new code matches old code
if(FALSE){
  library("lavaan")
  library("dr")
  
  HS.model <- 'visual  =~ x1 + x2 + x3
                   textual =~ x4 + x5 + x6
                   speed   =~ x7 + x8 + x9 '
     
  fit <- cfa(HS.model, data=HolzingerSwineford1939)
  fit2 <- cfa(HS.model, data=HolzingerSwineford1939, group="school")

  source("vuong_source.R")
  tmp1 <- calcAB(fit)
  source("vs_old.R")
  tmp2 <- calcAB(fit)
  all.equal(tmp1$A, tmp2$A)

  tmp2 <- get.ll(fit, HolzingerSwineford1939[,7:15])
  source("vuong_source.R")
  tmp1 <- get.ll(fit)
  ## This differs from the old code because the old code incorrectly took means
  ## to be zero if meanstructure==FALSE.  Instead, we should set the means
  ## at the sample estimates.
  all.equal(tmp1, tmp2)
}
