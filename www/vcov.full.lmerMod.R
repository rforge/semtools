############################################################################
## vcov.full.lmerMod: vcov for objects of class lmerMod, including random
##                    effect parameters.
## Authors: Ting Wang & Edgar C. Merkle
## Paper: Robust standard errors for linear mixed effects models estimated
##        via lme4
## License: GPL (≥ 2)
############################################################################
vcov.full.lmerMod <- function(object, ...) {
  if (class(object@resp)[1] != "lmerResp") stop("estfun.lmerMod() only works for lmer() models.")
  if (object@devcomp$dims[10]!=0) stop("estfun.lmerMod() only works for ML estimation.")
  
  ## minus inverse of Block 1: fixed parameter variance covariance matrix.
  fixvar <- vcov(object)
  
  ## preparation for Block 4:
  ## get all elements by getME and exclude multiple random effect models.
  parts <- getME(object, "ALL")
  
  ## shortcut, y-X\beta
  yXbe <- parts$y-parts$X %*% parts$beta
  ## total variance
  V <- (parts$Z %*% parts$Lambda %*% parts$Lambdat %*% parts$Zt + 
        diag(1, nrow(parts$X))) * (parts$sigma)^2
  
  ## length = (variance covariance parameter in G + residual variance.)
  uluti <- length(parts$theta)
  devV <- vector("list", (uluti+1))
  
  ## get d(G)/d(sigma), faciliated by d(Lambda). 
  devLambda <- vector("list", uluti)
  score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
  for (i in 1:uluti){
    devLambda[[i]] <- parts$Lambda
    devLambda[[i]][which(devLambda[[i]] != parts$theta[i])] <- 0
    devLambda[[i]][which(parts$Lambda == parts$theta[i])] <- 1
    devLambda[[i]] <- forceSymmetric(devLambda[[i]], uplo="L")
    devV[[i]] <- (parts$Z %*% devLambda[[i]]  %*% parts$Zt)
  }
  devV[[(uluti+1)]] <- diag(1,nrow(parts$X))
  
  ## Block 4: sigma's second derivative.
  ranvar <- matrix(NA, nrow = (uluti+1), ncol = (uluti+1))
  ## combinations (allow repeat) to indicate entries
  entries <- t(cbind(combn((uluti+1),2), matrix(rep(1:(uluti+1), each=2), 2, (uluti+1))))
  entries <- entries[order(entries[,1], entries[,2]),]
  
  for (i in 1:nrow(entries)){
    ranvar[lower.tri(ranvar, diag=TRUE)][i] <- as.numeric(((1/2) * sum(diag(solve(V) %*% 
                                               devV[[entries[i,1]]] %*% 
                                               solve(V) %*% devV[[entries[i,2]]]))))
  }
  ranvar <- solve(forceSymmetric(ranvar, uplo="L"))
  
  ## block 2 and block 3: second derivative of sigma and beta.
  varcov_beta <- matrix(NA, length(devV), length(parts$beta))
  for (j in 1 : (length(devV))){
    varcov_beta[j,] <- as.vector(-t(parts$X) %*%
                                 (solve(V) %*% (devV[[j]]) %*%
                                  solve(V)) %*% yXbe)
  }
  
  ## Organize full_varcov
  full_varcov <- -solve(rbind(cbind(-solve(fixvar), t(varcov_beta)),
                              cbind(varcov_beta, -solve(ranvar))))
  
  colnames(full_varcov) <- c(rownames(summary(object)$coefficients),
                             paste("cov",
                                   names(parts$theta),
                                   sep="_"), "residual")
  full_varcov
}
