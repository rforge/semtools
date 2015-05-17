## Calculate logintlik via laplace method. Bayes factor is based on
## the difference between two logintlik.

logintlik <- function(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("beta1[1]",p))],
                      data = c(choicesurerisk, tlist(tbeta1_1=0),
                               plist),
                      priorphi = "dunif",
                      n.iter = 5000){
  ## require packages  
  stopifnot(require("MCMCpack"))
  stopifnot(require("Matrix"))
  stopifnot(require("mvtnorm"))
  
  m1 <- jags.model(model, data, n.chains=3)
  update(m1, n.iter)
  m1res <- coda.samples(m1, paraest, n.iter = n.iter)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## all possible parameter names in bayeslv.Rnw. 
  paranam = c("lambda1", "lambda2", "lambda3","b0",
              "int1", "beta0", "beta1[1]","beta1[2]",
               "beta2[1]", "beta2[2]", "beta3[1]","beta3[2]",
              "beta4","beta5","beta6", "nu1", "nu2",
              "nu3","nu4", "Inv_cov1", "phi", 
              "Inv_sig_e1", "Inv_sig_ee", "Inv_sig_f1")
  pname <- names(m1med)
  rnums <- sapply(paranam, grep, pname, fixed=TRUE)
  
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and
  ## fixed parameters
  if (priorphi == "dwish") {
      m1cov <- m1cov[-c(rnums$Inv_cov1[3], 
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)]),
                    -c(rnums$Inv_cov1[3],
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                        m1med[rnums$lambda2]==1)])]
  } 
  if (priorphi == "gamma") {
      m1cov <- m1cov
  }
  if (priorphi == "dunif") {
     m1cov <- m1cov[-c(
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)]),
                    -c(
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)])]
  }

  
  ## log(f(theta*))
  logfthet <- ifelse(model!="2fsemt", sum(dnorm(m1med[
              rnums$lambda1[which(m1med[rnums$lambda1]!=0 &
                                  m1med[rnums$lambda1]!=1)]], 0,
              sqrt(1/0.001),log=TRUE)), sum(dnorm(m1med[
              rnums$lambda1[which(m1med[rnums$lambda1]!=0 &
                                  m1med[rnums$lambda1]!=1)]], 0,
              sqrt(1/as.numeric(plist[1])),log=TRUE))/0.5) + 
              ifelse(model!="2fsemt", sum(dnorm(m1med[
              rnums$lambda1[which(m1med[rnums$lambda2]!=0 &
                                  m1med[rnums$lambda2]!=1)]], 0,
              sqrt(1/0.001),log=TRUE)), sum(dnorm(m1med[
              rnums$lambda2[which(m1med[rnums$lambda2]!=0 &
                                  m1med[rnums$lambda2]!=1)]], 0,
              sqrt(1/as.numeric(plist[2])),log=TRUE))/0.5) + 
              sum(dnorm(m1med[rnums$lambda3], 0,
                        sqrt(1/as.numeric(plist[3])), log=TRUE)/0.5) +
              sum(dnorm(m1med[c(rnums$int1,      
              rnums$beta0, rnums$b0)], 0, sqrt(1/0.001), log=TRUE)) +
              sum(dnorm(m1med[c(rnums$'beta1[1]', rnums$'beta1[2]')], 0,
                        sqrt(1/as.numeric(plist[4])), log=TRUE)) +
              sum(dnorm(m1med[c(rnums$'beta2[1]', rnums$'beta2[2]')], 0,
                        sqrt(1/as.numeric(plist[5])), log=TRUE)) +
              sum(dnorm(m1med[c(rnums$'beta3[1]', rnums$'beta3[2]')], 0,
                        sqrt(1/as.numeric(plist[6])), log=TRUE)) +
              sum(dnorm(m1med[rnums$beta4], 0, sqrt(1/as.numeric(plist[7])), log=TRUE)) +
              sum(dnorm(m1med[rnums$beta5], 0, sqrt(1/as.numeric(plist[8])), log=TRUE)) +
              sum(dnorm(m1med[rnums$beta6], 0, sqrt(1/as.numeric(plist[9])), log=TRUE)) +
              sum(dnorm(m1med[rnums$nu1], 0, sqrt(1/as.numeric(plist[10])), log=TRUE)) +
              sum(dnorm(m1med[rnums$nu2], 0, sqrt(1/as.numeric(plist[11])), log=TRUE)) +
              sum(dnorm(m1med[rnums$nu3], 0, sqrt(1/as.numeric(plist[12])), log=TRUE)) +
              sum(dnorm(m1med[rnums$nu4], 0, sqrt(1/as.numeric(plist[13])), log=TRUE)) +
              ifelse(is.na(rnums$Inv_cov1[3]),
              0, log(dwish(matrix(c(m1med[rnums$Inv_cov1]), 2, 2),2, 
                  matrix(c(1,0,0,1),2,2)))) +
              sum(dunif(m1med[rnums$phi], -1, 1, log=TRUE)) +
              sum(dgamma(m1med[rnums$Inv_sig_f1], .1, .1, log=TRUE)) + 
              sum(dgamma(m1med[rnums$Inv_sig_e1], .1, .1,log=TRUE)) +
              sum(dgamma(m1med[rnums$Inv_sig_ee], .1, .1,log=TRUE)) 

  ## log-conditional likelihood, with thetas integrated out
    lcl <- 0
    if (priorphi == "dwish"){
      lambda <- matrix(c(m1med[rnums$lambda1],
                          m1med[rnums$lambda2]), 10, 2)
      phi <- solve(matrix(c(m1med[rnums$Inv_cov1]), 2, 2))
      psi <- diag(1/m1med[rnums$Inv_sig_e1])
      covmat <- lambda %*% phi %*% t(lambda) + psi
    }
    if (priorphi == "gamma"){
      lambda <- matrix(c(m1med[rnums$lambda1]), 10, 1)
      phi <- 1/m1med[rnums$Inv_sig_f1]
      psi <- diag(1/m1med[rnums$Inv_sig_e1])
      covmat <- lambda %*% phi %*% t(lambda) + psi
    }
    if (priorphi == "dunif"){
      psix <- diag(1/m1med[rnums$Inv_sig_ee])
      lambda <- matrix(c(m1med[rnums$lambda1],
                         m1med[rnums$lambda2]), 10, 2)
      psiy <- diag(1/m1med[rnums$Inv_sig_e1])
      phi <- matrix(c(1, m1med[rnums$phi], m1med[rnums$phi], 1),
                    2, 2)
      Mu <- matrix(NA, data$N, 2)
      numain <- c(ifelse(is.na(rnums$nu1[1]), 0,
                      m1med[rnums$nu1]),
                  ifelse(is.na(rnums$nu2[1]), 0,
                      m1med[rnums$nu2]))
      nuint <- matrix(c(ifelse(is.na(rnums$nu3[1]), 0,
                      m1med[rnums$nu3])*(data$Num - data$Numc), 
               ifelse(is.na(m1med[rnums$nu4[1]]), 0,
                      m1med[rnums$nu4])*
                  (data$Num - data$Numc)), nrow=data$N, ncol=2)
    }
  
    for(i in 1:data$N){
      if (is.na(m1med[rnums$Inv_sig_ee][1])){
        mnvec <- m1med[rnums$int1]
        lcl <- lcl + dmvnorm(data$y[i,], mnvec, covmat, log=TRUE)
      } else{
           Mu[i,1] <- (ifelse(is.na(m1med[rnums$'beta1[1]'[1]]), 0,
                       m1med[rnums$'beta1[1]'])*data$Frame[i]) + 
                (ifelse(is.na(m1med[rnums$'beta2[1]'[1]]), 0,
                       m1med[rnums$'beta2[1]'])*
                           (data$Num[i] - data$Numc)) + 
                      (ifelse(is.na(m1med[rnums$'beta3[1]'[1]]), 0,
                         m1med[rnums$'beta3[1]'])*data$Frame[i]*
                                 (data$Num[i] - data$Numc))
           Mu[i,2] <- (ifelse(is.na(m1med[rnums$'beta1[2]'[1]]), 0,
                       m1med[rnums$'beta1[2]'])*data$Frame[i]) + 
                (ifelse(is.na(m1med[rnums$'beta2[2]'[1]]), 0,
                       m1med[rnums$'beta2[2]'])*
                           (data$Num[i] - data$Numc)) + 
                      (ifelse(is.na(m1med[rnums$'beta3[2]'[1]]), 0,
                         m1med[rnums$'beta3[2]'])*data$Frame[i]*
                                 (data$Num[i] - data$Numc))
           
    cons <-  (ifelse(is.na(m1med[rnums$beta4[1]]), 0,
                             m1med[rnums$beta4])*data$Frame[i]) + 
                (ifelse(is.na(m1med[rnums$beta5[1]]), 0,
                     m1med[rnums$beta5])*(data$Num[i]-data$Numc)) + 
             (ifelse(is.na(m1med[rnums$beta6[1]]), 0,
          m1med[rnums$beta6])*data$Frame[i]*(data$Num[i]-data$Numc))
    xcovmat <- (m1med[rnums$lambda3] * (1 + numain%*%phi%*%numain +
                nuint[i,]%*%phi%*%nuint[i,])) %*%
               t(m1med[rnums$lambda3]) + psix
    xmnvec <- m1med[rnums$b0]+ m1med[rnums$lambda3] * 
             (numain%*%Mu[i,] + 
              nuint[i,]%*%Mu[i,] + cons)
    ycovmat <- lambda %*% phi %*% t(lambda) + psiy
    ymnvec <- m1med[rnums$beta0] + lambda%*%Mu[i,]
    mnvec <- c(xmnvec, ymnvec)
    covmat <- bdiag(xcovmat,ycovmat)
    covmat[(nrow(xcovmat)+1):(nrow(xcovmat)+nrow(ycovmat)), 
           1:nrow(xcovmat)] = 
               (lambda %*% c(numain%*%phi+nuint[i,]%*%phi)) %*% 
               t(m1med[rnums$lambda3])
    covmat[1:nrow(xcovmat), (nrow(xcovmat)+1):(nrow(xcovmat)+
    nrow(ycovmat))] = t(covmat[(nrow(xcovmat)+1):(nrow(xcovmat)+
                        nrow(ycovmat)),1:nrow(xcovmat)])
    
    lcl <- lcl + dmvnorm(cbind(data$X[i,],data$y[i,]), mnvec, 
                         as.matrix(covmat), log=TRUE)
     }
  }

  ## Eq (7) of lewis + raftery
  logintlik <- (nrow(m1cov)/2)*log(2*pi) + 
                .5*log(det(m1cov)) + logfthet + lcl
  logintlik
}

if(FALSE){
  source("bfcalc.R")
  source("jagsmodels.R")

}

#############################################################
### Calculate Bayes factor via path sampling
# bfcalc <- function(model, data, S=10, fname)
#  {
#   if(file.exists(fname)){
#     load(fname)
#   } else {
#    ## Arguments:
#    ## model (character): file name containing jags model
#    ## data (list): data used for jags model
#    ## S (numeric): number of t values
#    ## fname (character): file name to save ubar results
#    ## n.iter (numeric): number of iteration
#    ## burnin (numeric): number of burnin
#    mod3 <- vector("list", S)
#    mod3res <- vector("list", S)
#    ubar <- rep(NA, S)
#    data <- c(data, list(t = NA))
#    do.parallel <- require("parallel")
#    if (do.parallel){
#      ubar <- mclapply(1:(S+1), function(s){
#        data$t <- (s-1)*(1/S)
#        mod3[[s]] <- jags.model(file = model, data = data,
#                                n.chains = 3)
#        update(mod3[[s]])
#        ## iteration and burnin and (thin) can be changed.
#        mod3res[[s]] <- coda.samples(mod3[[s]], "u", n.iter=5000)
#        summary(mod3res[[s]])$statistics[1]},
#        mc.cores = round(.75*detectCores()),
#        mc.preschedule = FALSE)
#      ubar <- unlist(ubar)
#    } else {
#       for (s in 1:(S+1)){
#         data$t <- (s-1)*(1/S)
#         mod3[[s]] <- jags.model(file = model, data = data,
#                               n.chains = 3)
#         ## iteration and burnin and (thin) can be changed.
#         update(mod3[[s]])
#         mod3res[[s]] <- coda.samples(mod3[[s]], "u", n.iter=5000)
#         ubar[s] <- summary(mod3res[[s]])$statistics[1]
#       }
#    }
#    save(ubar, file=fname)
#  }
#    ## logBF
#    logBF=0
#    for(s in 1:S){
#      logBF=logBF+(ubar[s+1]+ubar[s])*(1/S)/2
#    }
#    logBF
#}

 
