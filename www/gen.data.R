"gen.data" <- function(sampsize,ses)
  {
    # Generates data from a factor analysis model that violates
    # measurement invariance.  Also generates a continuous auxiliary
    # variable related to the invariance (loosely called "age").
    
    # ses is number of SEs between group 1's parameters and group 2's
    # parameters.
    
    require("mvtnorm")

    # Could eventually manipulate this so break point is not directly
    # in the middle
    num.y <- round(sampsize/2,0)
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
    
    list(datmat=datmat, age=age)
  }
