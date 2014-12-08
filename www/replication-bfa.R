
###############################################
## Preliminary packages, data, files
library("knitr")
## global knitr options
opts_chunk$set(echo=FALSE, fig.path='figure/vuong-',fig.align='center',fig.show='hold',size='footnotesize', cache.path="cache/", warning=FALSE, message = FALSE)

## packages
library("rjags", quietly=TRUE)
library("mvtnorm")
library("MCMCpack")

## auxiliary code
source("jagsmodels.R")

## real-world data
load("peters08.rda")


###############################################
## Figure 1
par(mfrow=c(1,2))
shades <- paste("grey", c(100, 95, 85, 70, 60, 30))
ptshading <- shades[peters08$Lipkus - 5]

with(peters08[peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Positive frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
#abline(3.5,.625)

with(peters08[!peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Negative frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
#abline(3.5,.625)


###############################################
## Data creation
set.seed(1090)
surerisk <- list(y = peters08[,8:17], N = 108, J = 10)
sureriskreg <- list(y = peters08[,8:17], Frame = peters08$posframe,
                           Num = peters08$Lipkus, N = 108, J = 10)
choicesurerisk <- list(y = peters08[,8:17], N = 108, J = 10, K = 5,
                       X = peters08[,3:7],
                       Frame = peters08$posframe, 
                       Num = peters08$Lipkus,
                       Numc=mean(peters08$Lipkus))

###############################################
## Two-factor model estimation
if(file.exists("mod2res.rda")){
  load("mod2res.rda")
} else {
  mod2 <- jags.model("2f.jag", data = surerisk, n.chains = 3)
  update(mod2, 5000)
  mod2res <- coda.samples(mod2, c("lambda1", "lambda2","phi"), 
             n.iter=25000)
  save(mod2res, file="mod2res.rda")
}


###############################################
## Bayes factor for one- versus two-factor model
if(file.exists("logBF2f1f.rda")){
  load("logBF2f1f.rda")
} else {
  ## 2f model
  m1 <- jags.model("2f.jag", surerisk, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, c("lambda1", "lambda2[2:10]",
                            "int", "Inv_sig_e", "Inv_cov"), 
                      n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## remove fixed parameter
  m1med <- m1med[-30]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## fixed parameter
  m1cov <- m1cov[-c(3,30),-c(3,30)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[15:42], 0, sqrt(1/.001), log=TRUE)) 
              + log(dwish(matrix(c(m1med[1:4]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +
              sum(dgamma(m1med[5:14], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m1med[25:29], 0, m1med[30:33], 0, 
                     m1med[34:42]), 10, 2)
  phi <- solve(matrix(c(m1med[1:4]), 2, 2))
  psi <- diag(1/m1med[5:14])
  covmat <- lambda %*% phi %*% t(lambda) + psi

  for(i in 1:surerisk$N){
    mnvec <- m1med[15:24]
    lcl <- lcl + dmvnorm(surerisk$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(m1med)/2)*log(2*pi) + 
                .5*log(det(m1cov)) + logfthet + lcl
  
  ## Now repeat for 1f model
  m0 <- jags.model("1f.jag", surerisk, n.chains=3)
  update(m0, 1000)
  m0res <- coda.samples(m0, c("lambda","int", "Inv_sig_e", 
                                 "Inv_sig_f1"), n.iter=5000)
  ## posterior medians
  m0med <- summary(m0res)$quantiles[,3]
  ## posterior covariances
  m0cov <- cov(rbind(m0res[[1]], m0res[[2]], m0res[[3]]))
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m0med[12:31], 0, sqrt(1/.001), log=TRUE)) 
              + sum(dgamma(m0med[1:11], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m0med[22:31]), 10, 1)
  phi <- 1/m0med[11]
  psi <- diag(1/m0med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi

  for(i in 1:surerisk$N){
    mnvec <- m0med[12:21] 
    lcl <- lcl + dmvnorm(surerisk$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(m0med)/2)*log(2*pi) + 
                .5*log(det(m0cov)) + 
                logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBF2f1f <- logintlik1 - logintlik0
  save(logBF2f1f, file="logBF2f1f.rda")
 }  


###############################################
## Bayes factor for two-factor exploratory vs
## two-factor confirmatory model (not used in
## the paper)
if(file.exists("logBF2fc.rda")){
  load("logBF2fc.rda")
} else {
  ## 2f exploratory loading
  m1 <- jags.model("2f.jag", surerisk, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, c("lambda1", "lambda2[2:10]",
                            "int", "Inv_sig_e", "Inv_cov"), 
                        n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## remove fixed parameter
  m1med <- m1med[-30]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## fixed parameter
  m1cov <- m1cov[-c(3,30),-c(3,30)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[15:42], 0, sqrt(1/.001), log=TRUE)) 
              + log(dwish(matrix(c(m1med[1:4]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +
              sum(dgamma(m1med[5:14], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m1med[25:29], 0, m1med[30:33], 0, 
                     m1med[34:42]), 10, 2)
  phi <- solve(matrix(c(m1med[1:4]), 2, 2))
  psi <- diag(1/m1med[5:14])
  covmat <- lambda %*% phi %*% t(lambda) + psi

  for(i in 1:surerisk$N){
    mnvec <- m1med[15:24]
    lcl <- lcl + dmvnorm(surerisk$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(m1med)/2)*log(2*pi) + .5*log(det(m1cov)) + 
                logfthet + lcl
  
  ## Now repeat for 2f confirmatory loading
  m0 <- jags.model("2fre.jag", surerisk, n.chains=3)
  update(m0, 1000)
  m0res <- coda.samples(m0, c("lambda1[1:5]", "lambda2[6:10]",
                            "int", "Inv_sig_e", "Inv_cov"), 
                        n.iter=5000)
  ## posterior medians
  m0med <- summary(m0res)$quantiles[,3]
  ## posterior covariances
  m0cov <- cov(rbind(m0res[[1]], m0res[[2]], m0res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## fixed parameter
  m0cov <- m0cov[-c(3),-c(3)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m0med[15:34], 0, sqrt(1/.001), log=TRUE)) 
              + log(dwish(matrix(c(m0med[1:4]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +
              sum(dgamma(m0med[5:14], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m0med[25:29], rep(0,10),
                     m0med[30:34]), 10, 2)
  phi <- solve(matrix(c(m0med[1:4]), 2, 2))
  psi <- diag(1/m0med[5:14])
  covmat <- lambda %*% phi %*% t(lambda) + psi

  for(i in 1:surerisk$N){
    mnvec <- m0med[15:24]
    lcl <- lcl + dmvnorm(surerisk$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(m0med)/2)*log(2*pi) + .5*log(det(m0cov)) + 
                logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBF2fc <- logintlik1 - logintlik0
  save(logBF2fc, file="logBF2fc.rda")
 }  


###############################################
## Table 1
loadings <- summary(mod2res)$statistics[1:20,1]
loadmat <- t(matrix(round(loadings, 2), 10, 2))
loadtab <- apply(loadmat, 1, function(x) paste(x, collapse=" & "))
loadtab <- paste(c("Riskless & ", "Risky & "), 
                 loadtab, "\\\\", sep="")
writeLines(loadtab)


## Frequentist measurement invariance tests
## (shhh, don't tell anybody this is here)
## ## Check what lavaan has to say
## library(lavaan)
## m1 <- 'f1 =~ Sure1.1 + Sure1.2 + Sure1.3 + Sure1.4 + Sure1.5
##        f2 =~ Risky1.1 + Risky1.2 + Risky1.3 + Risky1.4 + Risky1.5'
## m1fit <- cfa(m1, peters08, group="posframe")
## m4fit <- cfa(m1, peters08, group="posframe",
##              group.equal=c("loadings","intercepts", "means"))
## library(semTools)
## measurementInvariance(m1, peters08, group="posframe")


###############################################
## Measurement invariance Bayes factor (equal
## loadings or not)
if(file.exists("logBFmi1.rda")){
  load("logBFmi1.rda")
} else {
  ## loading are different wrt frame
  srr <- c(surerisk, list(t=1))
  m1 <- jags.model("mi_loading.jag", srr, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, c("lambda1[2:5]", "lambda2[7:10]",
                              "delta1[2:5]", "delta2[7:10]",
                              "int", "int2", "Inv_sig_e", 
                              "Inv_sig_e2", "Inv_cov", "Inv_cov2"), 
                        n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] 
  m1cov <- m1cov[-c(3,7),-c(3,7)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[29:64], 0, sqrt(1/.001), log=TRUE)) + log(dwish(matrix(c(m1med[1:4]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +
               log(dwish(matrix(c(m1med[5:8]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +         
               sum(dgamma(m1med[9:28], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl1 <- 0
  lcl2 <- 0
  lambda1 <- matrix(c(1, m1med[57:60], rep(0,10), 1, m1med[61:64]),
                    10, 2)
  lambda2 <- matrix(c(1, m1med[57:60]+m1med[29:32], rep(0,10), 
                      1, m1med[61:64]+m1med[33:36]), 10, 2)
  phi1 <- solve(matrix(c(m1med[1:4]), 2, 2))
  phi2 <- solve(matrix(c(m1med[5:8]), 2, 2))
  psi1 <- diag(1/m1med[9:18])
  psi2 <- diag(1/m1med[19:28])
  covmat1 <- lambda1 %*% phi1 %*% t(lambda1) + psi1
  covmat2 <- lambda2 %*% phi2 %*% t(lambda2) + psi2
  
  for(i in 1:(surerisk$N/2-1)){
    mnvec <- m1med[37:46]
    lcl1 <- lcl1 + dmvnorm(surerisk$y[i,], mnvec, covmat1, log=TRUE)
  }
  for(i in (surerisk$N/2):surerisk$N){
    mnvec <- m1med[47:56]
    lcl2 <- lcl2 + dmvnorm(surerisk$y[i,], mnvec, covmat2, log=TRUE)
  }
  lcl <- lcl1 + lcl2
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(m1med)/2)*log(2*pi) + .5*log(det(m1cov)) + 
                logfthet + lcl
  
  ## Now repeat for equal loading model
  srr <- c(surerisk, list(t=0))
  m0 <- jags.model("mi_loading.jag", srr, n.chains=3)
  update(m0, 1000)
  m0res <- coda.samples(m0, c("lambda1[2:5]", "lambda2[7:10]",
                              "int", "int2", "Inv_sig_e", 
                              "Inv_sig_e2", "Inv_cov", "Inv_cov2"), 
                        n.iter=5000)
  ## posterior medians
  m0med <- summary(m0res)$quantiles[,3]
  ## posterior covariances
  m0cov <- cov(rbind(m0res[[1]], m0res[[2]], m0res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## fixed parameter
  m0cov <- m0cov[-c(3,7),-c(3,7)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m0med[29:56], 0, sqrt(1/.001), log=TRUE)) + log(dwish(matrix(c(m0med[1:4]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +
                 log(dwish(matrix(c(m0med[5:8]),2,2),2, 
                          matrix(c(1,0,0,1),2,2))) +                
              sum(dgamma(m0med[9:28], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl1 <- 0
  lcl2 <- 0
  lambda1 <- matrix(c(1, m0med[49:52], rep(0,10), 1, m0med[53:56]),
                    10, 2)
  phi1 <- solve(matrix(c(m0med[1:4]), 2, 2))
  phi2 <- solve(matrix(c(m0med[5:8]), 2, 2))
  psi1 <- diag(1/m0med[9:18])
  psi2 <- diag(1/m0med[19:28])
  covmat1 <- lambda1 %*% phi1 %*% t(lambda1) + psi1
  covmat2 <- lambda1 %*% phi2 %*% t(lambda1) + psi2
  
  for(i in 1:(surerisk$N/2-1)){
    mnvec <- m0med[29:38]
    lcl1 <- lcl1 + dmvnorm(surerisk$y[i,], mnvec, covmat1, log=TRUE)
  }
  for(i in (surerisk$N/2):surerisk$N){
    mnvec <- m0med[39:48]
    lcl2 <- lcl2 + dmvnorm(surerisk$y[i,], mnvec, covmat2, log=TRUE)
  }
  lcl <- lcl1 + lcl2
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(m0med)/2)*log(2*pi) + .5*log(det(m0cov)) + 
                logfthet + lcl
  
  ## log-bf of equal loading vs not equal
  logBFmi1 <- logintlik1 - logintlik0
  save(logBFmi1, file="logBFmi1.rda")
 }  


###############################################
## Measurement invariance Bayes factor (equal
## intercepts or not (equal loadings assumed)
if(file.exists("logBFmi2.rda")){
  load("logBFmi2.rda")
} else {
  srr <- c(surerisk, list(t=1))
  ## full model file mi_loadingintercept.jag.parameter changes. 
  m1 <- jags.model("mi_intercept.jag", srr, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, c("lambda1[2:5]", "lambda2[7:10]",
                              "int", "delta", "Inv_sig_e", 
                              "Inv_sig_e2", "Inv_cov", 
                              "Inv_cov2"), n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## Inv_cov2[1:2]=Inv_cov2[2:1]
  m1cov <- m1cov[-c(3,7),-c(3,7)]

  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[29:56], 0, sqrt(1/.001), log=TRUE)) + 
              log(dwish(matrix(m1med[1:4], 2, 2), 2, matrix(c(1,0,0,1), 2, 2))) +
              log(dwish(matrix(m1med[5:8], 2, 2), 2, matrix(c(1,0,0,1), 2, 2))) +
              sum(dgamma(m1med[9:28], .01, .01, log=TRUE))

  ## log-conditional likelihood, with thetas integrated out
  lcl1 <- 0
  lcl2 <- 0
  lambda <- matrix(c(1,m1med[49:52], rep(0,10), 1, m1med[53:56]), 10, 2)
  phi1 <- solve(matrix(m1med[1:4], 2, 2))
  phi2 <- solve(matrix(m1med[5:8], 2, 2))
  psi1 <- diag(1/m1med[9:18])
  psi2 <- diag(1/m1med[19:28])
  covmat1 <- lambda %*% phi1 %*% t(lambda) + psi1
  covmat2 <- lambda %*% phi2 %*% t(lambda) + psi2

  for(i in 1:(srr$N/2-1)){
    mnvec <- m1med[39:48]
    lcl1 <- lcl1 + dmvnorm(srr$y[i,], mnvec, covmat1, log=TRUE)
  }

  for(i in (srr$N/2):srr$N){
    mnvec <- m1med[39:48] + m1med[29:38]
    lcl2 <- lcl2 + dmvnorm(srr$y[i,], mnvec, covmat2, log=TRUE)
  }

  lcl <- lcl1 + lcl2

  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(m1med)/2)*log(2*pi) + .5*log(det(m1cov)) + 
                logfthet + lcl


  ## Now repeat for t=0
  ## Now repeat for equal loading+intercept model
  srr <- c(surerisk, list(t=0))
  m0 <- jags.model("mi_intercept.jag", srr, n.chains=3)
  update(m0, 1000)
  m0res <- coda.samples(m0, c("lambda1[2:5]", "lambda2[7:10]", "int", "Inv_sig_e", 
                              "Inv_sig_e2", "Inv_cov", "Inv_cov2"), 
                        n.iter=5000)
  ## posterior medians
  m0med <- summary(m0res)$quantiles[,3]
  ## posterior covariances
  m0cov <- cov(rbind(m0res[[1]], m0res[[2]], m0res[[3]]))
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and 
  ## Inv_cov2[1:2]=Inv_cov2[2:1]
  m0cov <- m0cov[-c(3,7),-c(3,7)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m0med[29:46], 0, sqrt(1/.001), log=TRUE)) + 
              log(dwish(matrix(m0med[1:4], 2, 2), 2, matrix(c(1,0,0,1), 2, 2))) +
              log(dwish(matrix(m0med[5:8], 2, 2), 2, matrix(c(1,0,0,1), 2, 2))) +
              sum(dgamma(m0med[9:28], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl1 <- 0
  lcl2 <- 0
  lambda1 <- matrix(c(1, m0med[39:42], rep(0,10), 1, m0med[43:46]), 10, 2)
  phi1 <- solve(matrix(c(m0med[1:4]), 2, 2))
  phi2 <- solve(matrix(c(m0med[5:8]), 2, 2))
  psi1 <- diag(1/m0med[9:18])
  psi2 <- diag(1/m0med[19:28])
  covmat1 <- lambda1 %*% phi1 %*% t(lambda1) + psi1
  covmat2 <- lambda1 %*% phi2 %*% t(lambda1) + psi2
  
  for(i in 1:(surerisk$N/2-1)){
    mnvec <- m0med[29:38]
    lcl1 <- lcl1 + dmvnorm(surerisk$y[i,], mnvec, covmat1, log=TRUE)
  }
  for(i in (surerisk$N/2):surerisk$N){
    mnvec <- m0med[29:38]
    lcl2 <- lcl2 + dmvnorm(surerisk$y[i,], mnvec, covmat2, log=TRUE)
  }
  lcl <- lcl1 + lcl2
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(m0med)/2)*log(2*pi) + .5*log(det(m0cov)) + 
                logfthet + lcl
  
  ## log-bf of equal loading+intercept vs not equal.
  logBFmi2 <- logintlik1 - logintlik0
  save(logBFmi2, file="logBFmi2.rda")
}


###############################################
## Two-factor model with covariates (no choice
## data included).  Not used in paper.
if(file.exists("summ2regres.rda")){
  load("summ2regres.rda")
} else {
    m2reg <- jags.model("2freg.jag", data = sureriskreg, 
                        n.chains = 3)
    update(m2reg,5000)
    m2regres <- coda.samples(m2reg, c("lambda1[1:5]",
                                      "lambda2[6:10]", 
                                      "beta0", "beta1", 
                                      "beta2", "beta3", "phi"),
                             n.iter=20000)
    summ2regres <- summary(m2regres)
    save(summ2regres, file="summ2regres.rda")
}
#summ2regres


###############################################
## Bayes factor for numeracy-by-frame interaction
## on riskless attraction
if(file.exists("logBFint.rda")){
  load("logBFint.rda")
} else {
  srr <- c(sureriskreg, list(t=1))
  m1 <- jags.model("linkinteraction.jag", srr, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, c("beta0", "lambda1[1:5]", 
                              "lambda2[6:10]", 
                              "beta1", "beta2", "beta3", 
                              "Inv_sig_e", "phi"), 
                        n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[11:34], 0, sqrt(1/.001), log=TRUE)) + 
              dunif(m1med[35], -1, 1, log=TRUE) +
              sum(dgamma(m1med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m1med[25:29], rep(0,10), m1med[30:34]), 10, 2)
  phi <- matrix(c(1, m1med[35], m1med[35], 1), 2, 2)
  psi <- diag(1/m1med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- m1med[11:20] + c(rep(m1med[21],5),rep(m1med[22],5)) * 
             srr$Frame[i] + 
             m1med[23]*srr$Num[i] + 
             rep(c(1,0), each=5)*m1med[24]*srr$Frame[i]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(m1med)/2)*log(2*pi) + .5*log(det(m1cov)) + 
                logfthet + lcl
  
  
  ## Now repeat for t=0
  srr$t <- 0
  m0 <- jags.model("linkinteraction.jag", srr, n.chains=3)
  update(m0, 1000)
  m0res <- coda.samples(m0, c("beta0", "lambda1[1:5]", 
                              "lambda2[6:10]", "beta1","beta2",
                              "Inv_sig_e", "phi"), n.iter=5000)
  ## posterior medians
  m0med <- summary(m0res)$quantiles[,3]
  ## posterior covariances
  m0cov <- cov(rbind(m0res[[1]], m0res[[2]], m0res[[3]]))
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m0med[11:33], 0, sqrt(1/.001), log=TRUE)) + 
              dunif(m0med[34], -1, 1, log=TRUE) +
              sum(dgamma(m0med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(m0med[24:28], rep(0,10), m0med[29:33]), 10, 2)
  phi <- matrix(c(1, m0med[34], m0med[34], 1), 2, 2)
  psi <- diag(1/m0med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- m0med[11:20] + c(rep(m0med[21],5), rep(m0med[22],5))*
             srr$Frame[i] + 
             m0med[23]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(m0med)/2)*log(2*pi) + .5*log(det(m0cov)) + 
                logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBFint <- logintlik1 - logintlik0
  save(logBFint, file="logBFint.rda") 
}


###############################################
## Bayes factor for frame effect on riskless
## attraction
if(file.exists("logBFsure1.rda")){
  load("logBFsure1.rda")
} else {
  srr <- c(sureriskreg, list(t=1))
  mbeta1 <- jags.model("linkFramesure1.jag", srr, n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1", 
                                     "lambda2", "beta1", "beta2",
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## remove fixed zeros
  beta1med <- beta1med[-(28:37)]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))
  beta1cov <- beta1cov[-(28:37), -(28:37)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[11:32], 0, sqrt(1/.001), 
              log=TRUE)) + 
              dunif(beta1med[33], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[23:27], rep(0,10), beta1med[28:32]), 
                   10, 2)
  phi <- matrix(c(1, beta1med[33], beta1med[33], 1), 2, 2)
  psi <- diag(1/beta1med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- beta1med[11:20] + c(rep(beta1med[21], 5), rep(0, 5))* 
             srr$Frame[i] + beta1med[22]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
                .5*log(det(beta1cov)) + logfthet + lcl
  
  
  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkFramesure1.jag", srr, n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1", 
                                     "lambda2", "beta2",
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## remove fixed zeros
  beta0med <- beta0med[-(27:36)]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))
  beta0cov <- beta0cov[-(27:36), -(27:36)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[11:31], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta0med[32], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[22:26], rep(0,10), beta0med[27:31]), 
                   10, 2)
  phi <- matrix(c(1, beta0med[32], beta0med[32], 1), 2, 2)
  psi <- diag(1/beta0med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- beta0med[11:20] + beta0med[21]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + 
                .5*log(det(beta0cov)) + logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBFsure1=logintlik1 - logintlik0
  save(logBFsure1, file="logBFsure1.rda")
  }


###############################################
## Bayes factor for frame effect on risky
## attraction
if(file.exists("logBFrisky1.rda")){
  load("logBFrisky1.rda")
} else {
  srr <- c(sureriskreg, list(t=1))
  mbeta1 <- jags.model("linkFramerisky1.jag", srr, n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1", 
                                     "lambda2", "beta1", "beta2",
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## remove fixed zeros
  beta1med <- beta1med[-(28:37)]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))
  beta1cov <- beta1cov[-(28:37), -(28:37)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[11:32], 0, sqrt(1/.001), 
              log=TRUE)) + 
              dunif(beta1med[33], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[23:27], rep(0,10), 
                     beta1med[28:32]), 10, 2)
  phi <- matrix(c(1, beta1med[33], beta1med[33], 1), 2, 2)
  psi <- diag(1/beta1med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- beta1med[11:20] + c(rep(0, 5), 
             rep(beta1med[21], 5)) * srr$Frame[i] + 
             beta1med[22]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
                .5*log(det(beta1cov)) + logfthet + lcl
  
  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkFramerisky1.jag", srr, n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1", 
                                     "lambda2", "beta2",
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## remove fixed zeros
  beta0med <- beta0med[-(27:36)]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))
  beta0cov <- beta0cov[-(27:36), -(27:36)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[11:31], 0, sqrt(1/.001), 
              log=TRUE)) + 
              dunif(beta0med[32], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:10], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[22:26], rep(0,10), beta0med[27:31]), 
                   10, 2)
  phi <- matrix(c(1, beta0med[32], beta0med[32], 1), 2, 2)
  psi <- diag(1/beta0med[1:10])
  covmat <- lambda %*% phi %*% t(lambda) + psi
  for(i in 1:srr$N){
    mnvec <- beta0med[11:20] + beta0med[21]*srr$Num[i]
    lcl <- lcl + dmvnorm(srr$y[i,], mnvec, covmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + 
                .5*log(det(beta0cov)) + logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBFrisky1 <- logintlik1 - logintlik0
  save(logBFrisky1, file="logBFrisky1.rda")
}


###############################################
## The main model, including two factors for
## attraction, choice data, numeracy, and so on
if(file.exists("summ2predres.rda")){
  load("summ2predres.rda")
} else {
    m2pred <- jags.model("2fpred.jag", data = choicesurerisk, 
                         n.chains = 3)
    update(m2pred,5000)
    m2predres <- coda.samples(m2pred, c("lambda1[1:5]", 
                                        "lambda2[6:10]", 
                                        "g1", "g2","g3","g4",
                                        "b0","b1","b2","b3",
                                        "beta0", "beta1", "beta2", 
                                        "beta3", "phi"), 
                              n.iter=20000)
    summ2predres <- summary(m2predres)
    save(summ2predres, file="summ2predres.rda")
}


## ----laplace-2fpred-interactionsure, echo=FALSE, eval=TRUE---------------
if(file.exists("logBFg4.rda")){
  load("logBFg4.rda")
} else {
  srr <- c(choicesurerisk, list(t=1))
  mbeta1 <- jags.model("linkinteractionchoicesure.jag", srr, 
                     n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g3","g4",  "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))

  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[16:46], 0, sqrt(1/.001), 
              log=TRUE)) + 
              dunif(beta1med[47], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:15], .01, .01, log=TRUE))

  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[37:41], rep(0,10), beta1med[42:46]), 
                   10, 2)
  phi <- matrix(c(1, beta1med[47], beta1med[47], 1), 2, 2)
  psi <- diag(1/beta1med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
           
  for(i in 1:srr$N){
    xcovmat <- beta1med[33]^2*covmat[6:10,6:10] + 
               beta1med[34]^2*covmat[1:5,1:5] - 
               beta1med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta1med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + 
               (beta1med[36]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] - 
               beta1med[35]*(srr$Num[i]-srr$Numc)*beta1med[36]*
               (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
               covmat[1:5,6:10]+ psi
  
    mnvec <- beta1med[16:20] + rep(beta1med[21],5)*srr$Frame[i] + 
             beta1med[22]*(srr$Num[i]-srr$Numc) + beta1med[33]*
             beta1med[28:32] + beta1med[34]*
             beta1med[23:27] + 
             beta1med[35]*beta1med[28:32]*
             (srr$Num[i]-srr$Numc) + 
             beta1med[36]*beta1med[23:27]* 
             (srr$Num[i]-srr$Numc)
  
  lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
}

  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
              .5*log(det(beta1cov)) + logfthet + lcl


  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkinteractionchoicesure.jag", srr, 
                       n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g3", "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))

  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[16:45], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta0med[46], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:15], .01, .01, log=TRUE))

  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[36:40], rep(0,10), 
                     beta0med[41:45]), 10, 2)
  phi <- matrix(c(1, beta0med[46], beta0med[46], 1), 2, 2)
  psi <- diag(1/beta0med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
  for(i in 1:srr$N){
    xcovmat <- beta0med[33]^2*covmat[6:10,6:10] + 
               beta0med[34]^2*covmat[1:5,1:5] - 
               beta0med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta0med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + psi
  
    mnvec <- beta0med[16:20] + rep(beta0med[21],5)*srr$Frame[i] +
             beta0med[22]*(srr$Num[i]-srr$Numc) + beta0med[33]*
             beta0med[28:32] + beta0med[34]*
             beta0med[23:27] + 
             beta0med[35]*beta0med[28:32]*
             (srr$Num[i]-srr$Numc) 
    
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }

  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + 
               .5*log(det(beta0cov)) + logfthet + lcl

  ## log-bf of frame vs without frame
  logBFg4 <- logintlik1 - logintlik0
  save(logBFg4, file="logBFg4.rda")
}


###############################################
## Bayes factor for "risky attraction by numeracy"
## interaction on choice
if(file.exists("logBFg3.rda")){
  load("logBFg3.rda")
} else {
  srr <- c(choicesurerisk, list(t=1))
  mbeta1 <- jags.model("linkinteractionchoicerisky.jag", srr, 
                     n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g3","g4",  "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))

  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[16:46], 0, sqrt(1/.001), 
              log=TRUE)) + 
              dunif(beta1med[47], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:15], .01, .01, log=TRUE))

  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[37:41], rep(0,10), beta1med[42:46]), 
                   10, 2)
  phi <- matrix(c(1, beta1med[47], beta1med[47], 1), 2, 2)
  psi <- diag(1/beta1med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
           
  for(i in 1:srr$N){
    xcovmat <- beta1med[33]^2*covmat[6:10,6:10] + 
               beta1med[34]^2*covmat[1:5,1:5] - 
               beta1med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta1med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + 
               (beta1med[36]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] - 
               beta1med[35]*(srr$Num[i]-srr$Numc)*beta1med[36]*
               (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
               covmat[1:5,6:10]+ psi
  
    mnvec <- beta1med[16:20] +  rep(beta1med[21],5)*srr$Frame[i] + 
             beta1med[22]*(srr$Num[i]-srr$Numc) + beta1med[33]*
             beta1med[28:32] + beta1med[34]*
             beta1med[23:27] + 
             beta1med[35]*beta1med[28:32]*
             (srr$Num[i]-srr$Numc) + 
             beta1med[36]*beta1med[23:27]* 
             (srr$Num[i]-srr$Numc)
  
  lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
}

  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
              .5*log(det(beta1cov)) + logfthet + lcl


  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkinteractionchoicerisky.jag", srr, 
                       n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g4", "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))

  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[16:45], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta0med[46], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:15], .01, .01, log=TRUE))

  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[36:40], rep(0,10), 
                     beta0med[41:45]), 10, 2)
  phi <- matrix(c(1, beta0med[46], beta0med[46], 1), 2, 2)
  psi <- diag(1/beta0med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
  for(i in 1:srr$N){
    xcovmat <- beta0med[33]^2*covmat[6:10,6:10] + 
               beta0med[34]^2*covmat[1:5,1:5] - 
               beta0med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta0med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] + psi
  
    mnvec <- beta0med[16:20] + rep(beta0med[21],5)*srr$Frame[i] +
             beta0med[22]*(srr$Num[i]-srr$Numc) + beta0med[33]*
             beta0med[28:32] + beta0med[34]*
             beta0med[23:27] + 
             beta0med[35]*beta0med[23:27]*
             (srr$Num[i]-srr$Numc) 
    
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }

  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + 
               .5*log(det(beta0cov)) + logfthet + lcl

  ## log-bf of frame vs without frame
  logBFg3 <- logintlik1 - logintlik0
  save(logBFg3, file="logBFg3.rda")
}


###############################################
## Bayes factor for the effect of frame on choice
if(file.exists("logBFfrchoice1.rda")){
  load("logBFfrchoice1.rda")
} else {
## Bayes factor based on linkFramesure.txt, using laplace
  srr <- c(choicesurerisk, list(t=1))
  mbeta1 <- jags.model("linkFramechoice1.jag", srr, n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1", "lambda2", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g3","g4",  "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## remove fixed zeros
  beta1med <- beta1med[-(42:51)]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))
  beta1cov <- beta1cov[-(42:51), -(42:51)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[16:46], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta1med[47], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:15], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[37:41], rep(0,10), 
                     beta1med[42:46]), 10, 2)
  phi <- matrix(c(1, beta1med[47], beta1med[47], 1), 2, 2)
  psi <- diag(1/beta1med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
             
  for(i in 1:srr$N){
    xcovmat <- beta1med[33]^2*covmat[6:10,6:10] + 
               beta1med[34]^2*covmat[1:5,1:5] - 
               beta1med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta1med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + 
               (beta1med[36]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] - 
               beta1med[35]*(srr$Num[i]-srr$Numc)*beta1med[36]*
               (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
               covmat[1:5,6:10]+ psi
    
    mnvec <- beta1med[16:20] +  beta1med[21]*srr$Frame[i] + 
             beta1med[22]*(srr$Num[i]-srr$Numc) + beta1med[33]*
             beta1med[28:32] + beta1med[34]*
             beta1med[23:27] + 
             beta1med[35]*beta1med[28:32]*
             (srr$Num[i]-srr$Numc) + 
             beta1med[36]*beta1med[23:27]* 
             (srr$Num[i]-srr$Numc)
    
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
                .5*log(det(beta1cov)) + logfthet + lcl
  
  
  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkFramechoice1.jag", srr, n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1", "lambda2", 
                                     "b0", "b2", "g1", "g2", 
                                     "g3","g4",  "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## remove fixed zeros
  beta0med <- beta0med[-(41:50)]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))
  beta0cov <- beta0cov[-(41:50), -(41:50)]
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[16:45], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta0med[46], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:15], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[36:40], rep(0,10), 
                     beta0med[41:45]), 10, 2)
  phi <- matrix(c(1, beta0med[46], beta0med[46], 1), 2, 2)
  psi <- diag(1/beta0med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
  for(i in 1:srr$N){
     xcovmat <- beta0med[32]^2*covmat[6:10,6:10] + 
               beta0med[33]^2*covmat[1:5,1:5] - 
               beta0med[32]*beta1med[33]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta0med[34]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + 
               (beta0med[35]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] - 
               beta0med[34]*(srr$Num[i]-srr$Numc)*beta0med[35]*
               (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
               covmat[1:5,6:10]+ psi
    
    mnvec <- beta0med[16:20]  + 
             beta0med[21]*(srr$Num[i]-srr$Numc) + beta0med[32]*
             beta0med[27:31] + beta0med[33]*
             beta0med[22:26] + 
             beta0med[34]*beta0med[27:31]*
             (srr$Num[i]-srr$Numc) + 
             beta0med[35]*beta0med[22:26]* 
             (srr$Num[i]-srr$Numc)
      
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + .5*log(det(beta0cov)) + logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBFfrchoice1 <- logintlik1 - logintlik0
  save(logBFfrchoice1, file="logBFfrchoice1.rda")
}


###############################################
## Bayes factor for the frame-by-numeracy
## interaction on choice.
if(file.exists("logBFb3choice.rda")){
  load("logBFb3choice.rda")
} else {
## Bayes factor based on linkFramesure.txt, using laplace
  srr <- c(choicesurerisk, list(t=1))
  mbeta1 <- jags.model("linkb3choice.jag", srr, n.chains=3)
  update(mbeta1, 1000)
  beta1res <- coda.samples(mbeta1, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "b3", "g1", 
                                     "g2", "g3","g4",  
                                     "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta1med <- summary(beta1res)$quantiles[,3]
  ## posterior covariances
  beta1cov <- cov(rbind(beta1res[[1]], beta1res[[2]], 
                        beta1res[[3]]))
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta1med[16:47], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta1med[48], -1, 1, log=TRUE) +
              sum(dgamma(beta1med[1:15], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta1med[38:42], rep(0,10), 
                     beta1med[43:47]), 10, 2)
  phi <- matrix(c(1, beta1med[48], beta1med[48], 1), 2, 2)
  psi <- diag(1/beta1med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
             
  for(i in 1:srr$N){
    xcovmat <- beta1med[34]^2*covmat[6:10,6:10] + 
               beta1med[35]^2*covmat[1:5,1:5] - 
               beta1med[34]*beta1med[35]*covmat[1:5,6:10] %*% 
               covmat[6:10,1:5] +
               (beta1med[36]*(srr$Num[i]-srr$Numc))^2*
               covmat[6:10,6:10] + 
               (beta1med[37]*(srr$Num[i]-srr$Numc))^2*
               covmat[1:5,1:5] - 
               beta1med[36]*(srr$Num[i]-srr$Numc)*beta1med[37]*
               (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
               covmat[1:5,6:10]+ psi
    
    mnvec <- beta1med[16:20] + rep(beta1med[21],5)*srr$Frame[i] + 
             beta1med[22]*(srr$Num[i]-srr$Numc) + 
             beta1med[23]*srr$Frame[i]*(srr$Num[i]-srr$Numc) +
             beta1med[34]*beta1med[29:33] + beta1med[35]*
             beta1med[24:28] + 
             beta1med[36]*beta1med[29:33]*
             (srr$Num[i]-srr$Numc) + 
             beta1med[37]*beta1med[24:28]* 
             (srr$Num[i]-srr$Numc)
    
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik1 <- (length(beta1med)/2)*log(2*pi) +
                .5*log(det(beta1cov)) + logfthet + lcl
  
  
  ## Now repeat for t=0
  srr$t <- 0
  mbeta0 <- jags.model("linkb3choice.jag", srr, n.chains=3)
  update(mbeta0, 1000)
  beta0res <- coda.samples(mbeta0, c("beta0", "lambda1[1:5]", 
                                     "lambda2[6:10]", 
                                     "b0", "b1", "b2", "g1", "g2", 
                                     "g3","g4", "Inv_sig_ee",  
                                     "Inv_sig_e", "phi"), 
                           n.iter=5000)
  ## posterior medians
  beta0med <- summary(beta0res)$quantiles[,3]
  ## posterior covariances
  beta0cov <- cov(rbind(beta0res[[1]], beta0res[[2]], 
                        beta0res[[3]]))
  
  ## log(f(theta*))
  logfthet <- sum(dnorm(beta0med[16:46], 0, sqrt(1/.001), 
                        log=TRUE)) + 
              dunif(beta0med[47], -1, 1, log=TRUE) +
              sum(dgamma(beta0med[1:15], .01, .01, log=TRUE))
  
  ## log-conditional likelihood, with thetas integrated out
  lcl <- 0
  lambda <- matrix(c(beta0med[37:41], rep(0,10), 
                     beta0med[42:46]), 10, 2)
  phi <- matrix(c(1, beta0med[47], beta0med[47], 1), 2, 2)
  psi <- diag(1/beta0med[11:15])
  covmat <- lambda %*% phi %*% t(lambda) 
  for(i in 1:srr$N){
     xcovmat <- beta0med[33]^2*covmat[6:10,6:10] + 
                beta0med[34]^2*covmat[1:5,1:5] - 
                beta0med[33]*beta1med[34]*covmat[1:5,6:10] %*% 
                covmat[6:10,1:5] +
                (beta0med[35]*(srr$Num[i]-srr$Numc))^2*
                covmat[6:10,6:10] + 
                (beta0med[36]*(srr$Num[i]-srr$Numc))^2*
                covmat[1:5,1:5] - 
                beta0med[35]*(srr$Num[i]-srr$Numc)*beta0med[36]*
                (srr$Num[i]-srr$Numc)*covmat[6:10,1:5] %*% 
                covmat[1:5,6:10]+ psi
    
     mnvec <- beta0med[16:20] + rep(beta0med[21],5)*srr$Frame[i]+
              beta0med[22]*(srr$Num[i]-srr$Numc) + beta0med[33]*
              beta0med[28:32] + beta0med[34]*
              beta0med[23:27] + 
              beta0med[35]*beta0med[28:32]*
              (srr$Num[i]-srr$Numc) + 
              beta0med[36]*beta0med[23:27]* 
              (srr$Num[i]-srr$Numc)
      
    lcl <- lcl + dmvnorm(srr$X[i,], mnvec, xcovmat, log=TRUE)
  }
  
  ## Eq (7) of lewis + raftery
  logintlik0 <- (length(beta0med)/2)*log(2*pi) + 
                .5*log(det(beta0cov)) + logfthet + lcl
  
  ## log-bf of frame vs without frame
  logBFb3choice <- logintlik1 - logintlik0
  save(logBFb3choice, file="logBFb3choice.rda")
}


###############################################
## Collect posterior intervals for Table 2;
## these are displayed with the above Bayes
## factors and F statistics from the original paper.
pname <- rownames(summ2predres$statistics)
rnums <- sapply(c("beta1", "b1", "beta3", "b3", "g3", "g4"), grep, 
                pname)
ints <- matrix(NA, 7, 2)
ints[1,] <- summ2predres$quantiles[rnums[[1]][1],c(1,5)]
ints[2,] <- summ2predres$quantiles[rnums[[1]][2],c(1,5)]
ints[3,] <- summ2predres$quantiles[rnums[[2]],c(1,5)]
ints[4,] <- summ2predres$quantiles[rnums[[3]],c(1,5)]
ints[5,] <- summ2predres$quantiles[rnums[[4]],c(1,5)]
ints[6,] <- summ2predres$quantiles[rnums[[5]],c(1,5)]
ints[7,] <- summ2predres$quantiles[rnums[[6]],c(1,5)]
ints <- apply(round(ints,2), 1, function(x) paste("(",paste(x, collapse=","),")",sep=""))
