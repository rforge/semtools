wdhfit <- function(data, silent = TRUE, suppressWarnings = TRUE, ...)
{
  ## data checks
  data <- as.data.frame(data)
  stopifnot(c("VERBREA","ABSREA","NUMER","gpa","grouping") %in% names(data))
  data <- data[, c("VERBREA","ABSREA","NUMER","gpa","grouping")]

  ## require package
  stopifnot(require("lavaan"))
  
  ## model from Step 5b of Wicherts et al., 2005:
  mod5b <- 'f1 =~ 1*VERBREA + label(rep("l2",4))*ABSREA + label(c("l1","l1","l1","l1expmin"))*NUMER 
                  f1 ~ label(c("a","faccontmin","facexpmaj","facexpmin"))*1 + c(0,NA,NA,NA)*1
                  f1 ~~ label(c("vmaj","vmin","vmaj","vmin"))*f1
                  NUMER ~~ label(c("e1","e1","e1","e1expmin"))*NUMER
                  ABSREA ~~ label(rep("e2",4))*ABSREA
                  VERBREA ~~ label(rep("e3",4))*VERBREA

                  NUMER ~ label(c("meanNA","meanNA","meanNA","meanNAexpmin"))*1
                  ABSREA ~ label(c("meanARcontmaj","meanARcontmin","meanAR","meanAR"))*1
                  VERBREA ~ label(rep("meanVR",4))*1'

  ## fit model
  rval <- lavaan(mod5b, data=data, meanstructure=TRUE, group="grouping")

  ## store (1) data, (2) fitted model object, (3) grouping variable,
  ## (4) focal group nubmers (group ordering is defined by lavaan and
  ##                          can be found in rval@Sample@group.label)
  rval <- list(
    data = data[,1:3],
    model = rval,
    groupvar = data$grouping,
    groupfoc = 4
  )
  class(rval) <- "wdhfit"

  return(rval)
}

wdhfit2 <- function(data, silent = TRUE, suppressWarnings = TRUE, ...)
{
  ## Fit the same model but with different identification constraints.
  ## data checks
  data <- as.data.frame(data)
  stopifnot(c("VERBREA","ABSREA","NUMER","gpa","grouping") %in% names(data))
  data <- data[, c("VERBREA","ABSREA","NUMER","gpa","grouping")]

  ## require package
  stopifnot(require("lavaan"))
  
  ## model from Step 5b of Wicherts et al., 2005:
  mod5b <- 'f1 =~ label(rep("l3",4))*VERBREA + 1*ABSREA + label(c("l1","l1","l1","l1expmin"))*NUMER 
                  f1 ~ label(c("a","faccontmin","facexpmaj","facexpmin"))*1 + c(0,NA,NA,NA)*1
                  f1 ~~ label(c("vmaj","vmin","vmaj","vmin"))*f1
                  NUMER ~~ label(c("e1","e1","e1","e1expmin"))*NUMER
                  ABSREA ~~ label(rep("e2",4))*ABSREA
                  VERBREA ~~ label(rep("e3",4))*VERBREA

                  NUMER ~ label(c("meanNA","meanNA","meanNA","meanNAexpmin"))*1
                  ABSREA ~ label(c("meanARcontmaj","meanARcontmin","meanAR","meanAR"))*1
                  VERBREA ~ label(rep("meanVR",4))*1'

  ## fit model
  rval <- lavaan(mod5b, data=data, meanstructure=TRUE, group="grouping")

  ## store (1) data, (2) fitted model object, (3) grouping variable,
  ## (4) focal group nubmers (group ordering is defined by lavaan and
  ##                          can be found in rval@Sample@group.label)
  rval <- list(
    data = data[,1:3],
    model = rval,
    groupvar = data$grouping,
    groupfoc = 4
  )
  class(rval) <- "wdhfit"

  return(rval)
}

logLik.wdhfit <- function(object, ...) {
  rval <- as.vector(fitMeasures(object$model, "logl"))

  structure(rval, df = length(coef(object)),
    nobs = nrow(object$data), class = "logLik")
}

vcov.wdhfit <- function(object, ...) {
  vc <- getMethod("vcov", "lavaan")(object$model, ...)
  class(vc) <- "matrix"
  vc
}

bread.wdhfit <- function(x, ...) vcov(x) * nrow(x$data)

info.wdhfit <- function(x, ...) solve(vcov(x) * nrow(x$data))

estfun.wdhfit <- function(x, ...)
{
  ## estimated coefficients
  cf <- coef(x)

  ## observed data
  X <- as.matrix(x$data)
  nvar <- ncol(X)

  ngrp <- ifelse(length(unique(x$groupvar)) > 1, x$model@Model@ngroups, 1)

  ## create numerical grouping variable
  ## needs groupvar object in x
  grpfac <- factor(x$groupvar, levels=unique(x$groupvar))
  levels(grpfac) <- 1:ngrp
  groupnum <- as.numeric(grpfac)
  groupnum.red <- groupnum[groupnum %in% x$groupfoc]
  nobs <- length(groupnum.red)

  ## Define matrix that will hold all scores
  ## Assumes means are estimated
  scores.H0 <- matrix(NA, nobs, x$model@Model@nx.free)

  for (i in x$groupfoc){
    ## fitted moments
    if(ngrp > 1){
      moments <- fitted(x$model)[[i]]
    }else{
      moments <- fitted(x$model)
    }
    Sigma.hat <- moments$cov
    Mu.hat <- moments$mean
    Sigma.inv <- solve(Sigma.hat)

    tmpX <- X[groupnum == i,]
    tmpnobs <- nrow(tmpX)

    ## Junk matrices for multiplication
    J <- matrix(1, 1, tmpnobs)
    J2 <- matrix(1, nvar, nvar); diag(J2) <- .5

    ## scores.H1 (H1 = saturated model)
    mean.diff <- t(t(tmpX) - Mu.hat %*% J)

    dx.Mu <- -1 * mean.diff %*% Sigma.inv

    dx.Sigma <- t(apply(mean.diff,1,function(x) lavaan:::vech(-J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat) %*% Sigma.inv))))

    scores.H1 <- cbind(dx.Mu, dx.Sigma)
    
    ## scores.H0
    Delta <- lavaan:::computeDelta(x$model@Model)[[i]]
    scores.H0[groupnum.red==i,] <- -scores.H1 %*% Delta
  }

  ## assign nice names and return
  rval <- scores.H0
  colnames(rval) <- names(cf)
  rownames(rval) <- rownames(X)[groupnum %in% x$groupfoc]
  return(rval)
}
