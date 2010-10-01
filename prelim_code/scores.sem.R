"scores.sem" <- function(sem.object,raw.data)
  {
    # Takes sem object as input, extracts scores (first derivatives
    # for each individual) using numerical derivatives and raw.data.
    # Code works but has not been fully tested.
    
    # The numerical derivatives are very slow, so that analytic
    # derivatives may be useful at some point.
    require(numDeriv)

    # Set up casewise likelihood function to be used in call to grad().
    # Does not currently work for models with intercepts
    # (Must worry about fixed.x argument for sem(); the means
    #  argument below would then go away):
    case.lik <- function(x, sem.object, ind.dat, means){
      # Get parameter locations in RAM matrices:
      param.locs <- sem.object$ram[sem.object$par.posn,]
      # Separate rows with heads==0 from rows with heads==1:
      A.locs <- param.locs[param.locs[,1]==1,]
      P.locs <- param.locs[param.locs[,1]==2,]

      # Fill in A and P matrices using tmp:
      mat.dim <- nrow(sem.object$A)
      A.tmp <- matrix(0,mat.dim,mat.dim)
      P.tmp <- matrix(0,mat.dim,mat.dim)
      for (i in 1:nrow(A.locs)){
        tmp.loc <- A.locs[i,2:3]
        A.tmp[tmp.loc[1],tmp.loc[2]] <- x[i]
      }
      for (i in 1:nrow(P.locs)){
        tmp.loc <- P.locs[i,2:3]
        P.tmp[tmp.loc[1],tmp.loc[2]] <- x[(nrow(A.locs) + i)]
      }
      
      # Model-implied S:
      S.imp <- with(sem.object, J %*% solve(diag(mat.dim) - A.tmp) %*%
                    P.tmp %*% t(solve(diag(mat.dim) - A.tmp)) %*% t(J))

      # Casewise likelihood:
      fval <- log(det(S.imp)) + t(ind.dat - means) %*% solve(S.imp) %*%
        (ind.dat - means)

      fval
    }

    # Use the function above to extract scores:
    scores <- matrix(NA,nrow(raw.data),length(sem.object$par.posn))
    means <- apply(raw.data,2,mean)
    for (j in 1:nrow(raw.data)){
      scores[j,] <- grad(case.lik, method="Richardson",
                         x=sem.object$coeff, sem.object=sem.object,
                         ind.dat=as.matrix(raw.data)[j,], means=means)
    }
    colnames(scores) <- rownames(sem.object$ram[sem.object$par.posn,])
    
    scores
  }
