"score.fa" <- function(res.list,derivs="a")
  {
    # Code to extract scores from a 2-factor model that
    # is fit by OpenMx.  res is a list containing
    # the model fit by mxRun() and the observed data matrix.

    # Possible options for derivs:
    #   a: Analytic derivatives (fast)
    #   n: Numerical derivatives (slow)

    derivs <- tolower(derivs)
    if (derivs=="n"){
      require("numDeriv")
    }

    res <- res.list$res
    dat <- res.list$dat

    # Results from OpenMx
    theta <- res@output$estimate

    # Get n and k from OpenMx output
    n <- nrow(dat)
    k <- length(theta)

    # Score matrix
    scores <- matrix(NA,n,k)
    
    switch(derivs, `n` = {
    
      # Define casewise likelihood function
      likefun <- function(x,ind.dat){
        mu <- x[14:19]
        lam <- matrix(c(x[1:3],rep(0,6),x[4:6]),6,2) # loadings
        phi <- matrix(c(1,x[13],x[13],1),2,2) # factor covs
        psi <- diag(x[7:12]) # error vars

        impl.cov <- lam %*% phi %*% t(lam) + psi

        # Added (-1/2) term:
        fval <- (-1/2)*(log(det(impl.cov)) +
                        t(ind.dat - mu) %*% solve(impl.cov) %*% (ind.dat - mu))
        fval
      }


      for (j in 1:n){
        scores[j,] <- grad(likefun, theta,
                           ind.dat=as.matrix(dat)[j,],
                           method="Richardson")
      }
    }, `a` = {
      mu <- theta[14:19]
      Lambda <- matrix(c(theta[1:3],rep(0,6),theta[4:6]),6,2) # loadings
      Phi <- matrix(c(1,theta[13],theta[13],1),2,2) # factor covs
      Psi <- matrix(0,6,6)
      diag(Psi) <- theta[7:12] # error vars

      Sigma <- Lambda %*% Phi %*% t(Lambda) + Psi
  
      Sig.inv <- solve(Sigma)
      for (i in 1:n){
        Xmu <- (dat[i,]-mu)
        SigiXmu <- Sig.inv %*% Xmu %*% t(Xmu) %*% Sig.inv

        dFdmu <- -as.vector(Sig.inv %*% Xmu)
        dFdlam <- as.vector((Sig.inv - SigiXmu) %*% Lambda %*% Phi)[c(1:3,10:12)]
        dFdphi <- (t(Lambda) %*% (Sig.inv - SigiXmu) %*% Lambda)[1,2]
        dFdpsi <- (1/2)*diag(Sig.inv - SigiXmu)
    
        scores[i,] <- -c(dFdlam, dFdpsi, dFdphi, dFdmu)
      }
    })

    scores
  }
           

"inf.fa" <- function(res.list, order.by, data)
  {
    # Code to extract information matrix of model parameters,
    # for strucchange.
    res <- res.list$res
    
    n <- nrow(res.list$dat)
    inf.mat <- (1/(2*n))*res@output$calculatedHessian

    inf.mat
  }
