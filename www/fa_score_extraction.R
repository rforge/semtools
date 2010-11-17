"score.fa" <- function(res.list)
  {
    # Code to extract scores from a 2-factor model that
    # is fit by OpenMx.  res is a list containing
    # the model fit by mxRun() and the observed data matrix.

    require(numDeriv)

    res <- res.list$res
    dat <- res.list$dat
    
    # Define casewise likelihood function
    likefun <- function(x,ind.dat){
      mu <- x[14:19]
      lam <- matrix(c(x[1:3],rep(0,6),x[4:6]),6,2) # loadings
      phi <- matrix(c(1,x[13],x[13],1),2,2) # factor covs
      psi <- diag(x[7:12]) # error vars

      impl.cov <- lam %*% phi %*% t(lam) + psi

      fval <- log(det(impl.cov)) + t(ind.dat - mu) %*% solve(impl.cov) %*%
                                   (ind.dat - mu)
      fval
    }

    # Results from OpenMx
    ests <- res@output$estimate
    # Need n from OpenMx output
    n <- nrow(dat)

    scores <- matrix(NA,n,length(ests))
    for (j in 1:n){
      scores[j,] <- grad(likefun, ests,
                         ind.dat=as.matrix(dat)[j,],
                         method="Richardson")
    }

    scores
  }


"inf.fa" <- function(res.list, order.by, data)
  {
    # Code to extract information matrix of model parameters,
    # for strucchange.
    res <- res.list$res
    
    n <- nrow(res.list$dat)
    inf.mat <- (2/n)*res@output$calculatedHessian

    inf.mat
  }
