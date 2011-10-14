mzfit <- function(data, engine = c("lavaan", "OpenMx"),
  silent = TRUE, suppressWarnings = TRUE, ...)
{
  ## data checks
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "y1", "y2", "y3") %in% names(data))
  data <- data[, c("x1", "x2", "x3", "y1", "y2", "y3")]

  ## engine for factor modeling
  engine <- match.arg(engine)

  switch(engine,
  
  "OpenMx" = {
    ## require package
    stopifnot(require("OpenMx"))

    ## model from Merkle & Zeileis
    rval <- mxModel("Factor Model", 
      type = "RAM",
      mxData(observed = data, type = "raw"),
      manifestVars = c("x1", "x2", "x3", "y1", "y2", "y3"), 
      latentVars = c("F1", "F2"),
      mxPath(from = c("x1", "x2", "x3", "y1", "y2", "y3"),
        arrows = 2, free = TRUE, values = c(1, 1, 1, 1, 1, 1),
        labels = c("e1", "e2", "e3", "e4", "e5", "e6"), lbound=0),
      mxPath(from = c("F1", "F2"),
        arrows = 2, all = TRUE, free = c(FALSE, TRUE, TRUE, FALSE),
        values = c(1, 0.5, 0.5, 1), labels = c("varF1", "cov", "cov", "varF2"),
        lbound=c(NA, 0, 0, NA)),
      mxPath(from = "F1", to = c("x1", "x2", "x3"),
        arrows = 1, free = TRUE, values = c(1, 1, 1),
        labels = c("l1", "l2", "l3"), lbound=0),
      mxPath(from = "F2", to = c("y1", "y2", "y3"),
        arrows = 1, free = TRUE, values = c(1, 1, 1),
        labels = c("l4", "l5", "l6"), lbound=0), 
      mxPath(from = "one",
        to = c("x1", "x2", "x3", "y1", "y2", "y3", "F1", "F2"),
        arrows = 1, free = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
        values = c(1, 1, 1, 1, 1, 1, 0, 0), labels = c("meanx1", "meanx2",
          "meanx3", "meany1", "meany2", "meany3", NA, NA),
        lbound = 0)
    )

    ## fit model
    rval <- mxRun(rval, silent = silent, suppressWarnings = suppressWarnings, ...)
  },
  
  "lavaan" = {
    ## require package
    stopifnot(require("lavaan"))

    ## model from Merkle & Zeileis
    rval <- 'verbal =~ x1 + x2 + x3; math =~ y1 + y2 + y3'

    ## fit model (set likelihood = "wishart" for N-1 correction as in OpenMx)
    rval <- sem(rval, data = data, meanstructure = TRUE, std.lv = TRUE)
  })

  ## store (1) fitted model object, (2) engine used,
  ## (3) match between OpenMx and lavaan coefficient labeling and ordering
  rval <- list(
    data = data,
    model = rval,
    engine = engine,
    names = cbind(c("l1", "l2", "l3", "l4", "l5", "l6", "e1", "e2", "e3", "e4", "e5", "e6",
        "cov", "meanx1", "meanx2", "meanx3", "meany1", "meany2", "meany3"),
      c("verbal=~x1", "verbal=~x2", "verbal=~x3", "math=~y1", "math=~y2", "math=~y3", "x1~~x1", "x2~~x2", "x3~~x3",
        "y1~~y1", "y2~~y2", "y3~~y3", "verbal~~math", "x1~1", "x2~1", "x3~1", "y1~1", "y2~1", "y3~1"))
  )
  class(rval) <- "mzfit"

  return(rval)
}

coef.mzfit <- function(object, ...) {
  switch(object$engine,  
  "OpenMx" = {
    n <- object$names
    cf <- object$model@output$estimate
    cf <- cf[n[,1]]
    names(cf) <- n[,2]
    cf
  },
  "lavaan" = {
    unclass(getMethod("coef", "lavaan")(object$model, ...))
  })
}

logLik.mzfit <- function(object, ...) {
  rval <- switch(object$engine,  
  "OpenMx" = {
    -0.5 * object$model@output$minimum
  },
  "lavaan" = {
    as.vector(fitMeasures(object$model, "logl"))
  })
  structure(rval, df = length(coef(object)),
    nobs = nrow(object$data), class = "logLik")
}

vcov.mzfit <- function(object, ...) {
  switch(object$engine,  
  "OpenMx" = {
    n <- object$names
    vc <- 2 * solve(object$model@output$calculatedHessian)
    vc <- vc[n[,1], n[,1]]
    colnames(vc) <- rownames(vc) <- n[,2]
    vc
  },
  "lavaan" = {
    vc <- getMethod("vcov", "lavaan")(object$model, ...)
    class(vc) <- "matrix"
    vc
  })  
}

bread.mzfit <- function(x, ...) vcov(x) * nrow(x$data)

info.mzfit <- function(x, ...) solve(vcov(x) * nrow(x$data))

estfun.mzfit <- function(x, ...)
{
  ## estimated coefficients
  cf <- coef(x)

  ## observed data
  X <- as.matrix(x$data)
  nvar <- ncol(X)
  nobs <- nrow(X)

  rval <- switch(x$engine,
  
  ## for lavaan we could in principle also handle
  ## other models:
  "lavaan" = {
    # fitted moments
    moments <- fitted(x$model)
    Sigma.hat <- moments$cov
    Mu.hat <- moments$mean
    Sigma.inv <- solve(Sigma.hat)

    # scores.H1 (H1 = saturated model)
    pstar <- nvar * (nvar + 1)/2 + nvar
    scores.H1 <- matrix(NA, nobs, pstar)

    # compute H1 scores per case
    for(i in 1:nrow(X)) {
      diff.i <- t(X[i,] - Mu.hat)
      dx.Mu <- -1 * t(diff.i %*% Sigma.inv)
      dx.Sigma <- (-1 * (Sigma.inv %*% (crossprod(diff.i) - Sigma.hat) %*% Sigma.inv))
      # correction for symmetry
      diag(dx.Sigma) <- diag(dx.Sigma)/2
      # in lavaan: first the means, then the covariances
      scores.H1[i,] <- c(dx.Mu, lavaan:::vech(dx.Sigma))
    }

    # scores.H0
    Delta <- lavaan:::computeDelta(x$model@Model)[[1]]
    scores.H0 <- scores.H1 %*% Delta
    -scores.H0
  },

  ## for OpenMx the particular Merkle & Zeileis model
  ## is hard-coded
  "OpenMx" = {

    # score matrix
    scores <- matrix(NA, nobs, length(cf))
  
    mu <- cf[14:19]
    Lambda <- matrix(c(cf[1:3], rep(0, 6), cf[4:6]), 6, 2) # loadings
    Phi <- matrix(c(1, cf[13], cf[13], 1), 2, 2)	   # factor covs
    Psi <- matrix(0, 6, 6)
    diag(Psi) <- cf[7:12]				   # error vars

    Sigma <- Lambda %*% Phi %*% t(Lambda) + Psi  
    Sig.inv <- solve(Sigma)
  
    for (i in 1:nobs){
      Xmu <- (X[i,] - mu)
      SigiXmu <- Sig.inv %*% Xmu %*% t(Xmu) %*% Sig.inv
      dFdmu <- -as.vector(Sig.inv %*% Xmu)
      dFdlam <- as.vector((Sig.inv - SigiXmu) %*% Lambda %*% Phi)[c(1:3, 10:12)]
      dFdphi <- (t(Lambda) %*% (Sig.inv - SigiXmu) %*% Lambda)[1,2]
      dFdpsi <- (1/2) * diag(Sig.inv - SigiXmu)  
      scores[i,] <- -c(dFdlam, dFdpsi, dFdphi, dFdmu)
    }

    scores
  })

  ## assign nice names and return
  colnames(rval) <- names(cf)
  rownames(rval) <- rownames(X)
  return(rval)
}
