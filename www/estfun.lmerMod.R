############################################################################
## estfun.lmerMod: Compute scores for objects of class lmerMod
## Authors: Ting Wang & Edgar C. Merkle
## Paper: Robust standard errors for linear mixed effects models estimated
##        via lme4
## License: GPL (≥ 2)
############################################################################
estfun.lmerMod <- function(object, level = 2,...) {
  if (class(object@resp)[1] != "lmerResp") stop("estfun.lmerMod() only works for lmer() models.")
  if (object@devcomp$dims[10]!=0) stop("estfun.lmerMod() only works for ML estimation.")
    
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  if (length(parts$l_i) > 1 & level == 2) stop("Multiple cluster variables detected. Supply 'level=1' argument to estfun.lmerMod().")
  
  ## prepare shortcuts
  yXbe <- parts$y - parts$X %*% parts$beta
  V <- (parts$Z %*% parts$Lambda %*% parts$Lambdat %*% parts$Zt + 
          diag(1, nrow(parts$X))) * (parts$sigma)^2
  
  ## adapt from Stroup book page 131, last eq,
  ## score for fixed effects parameter: score_beta=XR^{-1}(y-Zb-X\beta)
  score_beta <- t(t(parts$X) %*% solve(V)) * (yXbe)
  
  ## prepare for score of variance covariance parameters. 
  ## get d(G)/d(sigma), faciliated by d(Lambda).
  uluti <- length(parts$theta)
  devLambda <- vector("list", uluti)
  score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
  devV <- vector ("list", (uluti+1))
  
  for (i in 1:uluti){
    devLambda[[i]] <- parts$Lambda
    devLambda[[i]][which(devLambda[[i]] != parts$theta[i])] <- 0
    devLambda[[i]][which(parts$Lambda == parts$theta[i])] <- 1
    devLambda[[i]] <- forceSymmetric(devLambda[[i]], uplo="L")
    devV[[i]] <- (parts$Z %*% devLambda[[i]] %*% parts$Zt)
  }
  devV[[(uluti+1)]] <- diag(1,nrow(parts$X))
  
  ## score for variance covariance parameter
  score_varcov <- matrix(NA, nrow = nrow(parts$X), ncol = (uluti+1))
  for (j in 1:length(devV)){
    score_varcov[,j] <- as.vector(-(1/2) * diag(solve(V) %*% (devV[[j]])) +
        t((1/2) * t(parts$y-parts$X %*% parts$beta) %*%
          solve(V) %*% (devV[[j]]) %*% solve(V)) *
        (parts$y-parts$X %*% parts$beta))
  }
  
  ## Organize score matrix
  score <- cbind(as.matrix(score_beta), score_varcov)
  colnames(score) <- c(rownames(summary(object)$coefficients),
                       paste("cov",
                             names(parts$theta),
                             sep="_"), "residual")
  ## Clusterwise scores if level==2
  if (level == 2){
    index <- rep(1:parts$l_i, (parts$n/parts$l_i))
    index <- index[order(index)]
    score <- aggregate(x=score, by = list(index), FUN=sum)[,-1]
  }

  score
}
