mzfit <- function(data, silent = TRUE, suppressWarnings = TRUE, ...)
{
  ## data checks
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "y1", "y2", "y3") %in% names(data))
  data <- data[, c("x1", "x2", "x3", "y1", "y2", "y3")]

  ## need OpenMx for model fitting
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
  return(rval)
}

coef.MxModel <- function(object, ...) object@output$estimate

logLik.MxModel <- function(object, ...) {
  structure(-0.5 * object@output$minimum, df = length(object@output$estimate),
    nobs = object@data@numObs, class = "logLik")
}

vcov.MxModel <- function(object, ...)
  2 * solve(object@output$calculatedHessian)

bread.MxModel <- function(x, ...)
  2 * solve(x@output$calculatedHessian) * x@data@numObs

estfun.MxModel <- function(x, type = c("analytic", "numeric"), ...)
{
  ## numeric or analytic derivatives
  type <- match.arg(type)
  if(type == "numeric") stopifnot(require("numDeriv"))

  ## estimated coefficients
  theta <- coef(x)

  ## extract observed data (FIXME: for some reason reordered!)
  dat <- x@data@observed
  dat <- dat[order(as.numeric(rownames(dat))),]
  dat <- as.matrix(dat)
  n <- nrow(dat)

  # score matrix
  scores <- matrix(NA, n, length(theta))
  colnames(scores) <- names(theta)
  rownames(scores) <- rownames(dat)
  
  switch(type,
  
  "numeric" = {    
    # casewise likelihood function
    likefun <- function(x, ind.dat){
      mu <- x[14:19]
      lam <- matrix(c(x[1:3], rep(0,6), x[4:6]), 6, 2) # loadings
      phi <- matrix(c(1,x[13], x[13], 1), 2, 2)        # factor covs
      psi <- diag(x[7:12])			       # error vars

      impl.cov <- lam %*% phi %*% t(lam) + psi

      fval <- (-1/2) * (log(det(impl.cov)) + t(ind.dat - mu) %*% solve(impl.cov) %*% (ind.dat - mu))
      fval
    }

    for(j in 1:n) {
      scores[j,] <- grad(likefun, theta, ind.dat = as.matrix(dat)[j, ], method = "Richardson")
    }
  },
  
  "analytic" = {
    mu <- theta[14:19]
    Lambda <- matrix(c(theta[1:3], rep(0, 6), theta[4:6]), 6, 2) # loadings
    Phi <- matrix(c(1, theta[13], theta[13], 1), 2, 2)           # factor covs
    Psi <- matrix(0, 6, 6)
    diag(Psi) <- theta[7:12]                                     # error vars

    Sigma <- Lambda %*% Phi %*% t(Lambda) + Psi  
    Sig.inv <- solve(Sigma)
    
    for (i in 1:n){
      Xmu <- (dat[i,] - mu)
      SigiXmu <- Sig.inv %*% Xmu %*% t(Xmu) %*% Sig.inv
      dFdmu <- -as.vector(Sig.inv %*% Xmu)
      dFdlam <- as.vector((Sig.inv - SigiXmu) %*% Lambda %*% Phi)[c(1:3, 10:12)]
      dFdphi <- (t(Lambda) %*% (Sig.inv - SigiXmu) %*% Lambda)[1,2]
      dFdpsi <- (1/2) * diag(Sig.inv - SigiXmu)  
      scores[i,] <- -c(dFdlam, dFdpsi, dFdphi, dFdmu)
    }
  })

  return(scores)
}


dgp <- function(n = 200, diff = 3)
{
  # Generates data from a factor analysis model that violates
  # measurement invariance.  Also generates a continuous auxiliary
  # variable related to the invariance (loosely called "age").
  
  # ses is number of SEs between group 1's parameters and group 2's
  # parameters.
  sampsize <- n
  ses <- diff
  
  stopifnot(require("mvtnorm"))

  # Could eventually manipulate this so break point is not directly
  # in the middle
  num.y <- round(sampsize/2, 0)
  num.o <- sampsize - num.y

  age <- c(runif(num.y,13,16),runif(num.o,16,18))
  
  # Define parameter vectors/matrices:
  mu <- c(29.32,24.70,14.84,10.59,19.30,18.01)

  # Loadings for "young" individuals:
  lambda.y <- matrix(0,6,2)
  lambda.y[,1] <- c(4.92,2.96,5.96,0,0,0)
  lambda.y[,2] <- c(0,0,0,3.24,4.32,7.21)

  # Loadings for "old" individuals:
  lambda.o <- matrix(0,6,2)
  lambda.o[,1] <- lambda.y[,1] - (ses * c(.87,.6,1.32,rep(0,3)))
  lambda.o[,2] <- c(0,0,0,3.24,4.32,7.21)
  
  phi <- matrix(.48,2,2)
  diag(phi) <- 1

  psi <- matrix(0,6,6)
  diag(psi) <- c(26.77,13.01,30.93,3.17,8.82,22.5)

  # Initialize data matrix:
  datmat <- matrix(0,sampsize,6)
  colnames(datmat) <- c("x1", "x2", "x3", "y1", "y2", "y3")

  # Generate z and u vectors:
  z <- t(rmvnorm(sampsize,rep(0,2),phi))
  u <- t(rmvnorm(sampsize,rep(0,6),psi))

  # Based on z and u, generate data:
  for (i in 1:num.y){
    datmat[i,] <- mu + lambda.y%*%z[,i] + u[,i]
  }
  for (i in (num.y+1):sampsize){
    datmat[i,] <- mu + lambda.o%*%z[,i] + u[,i]
  }
  
  cbind(as.data.frame(datmat), age)
}
