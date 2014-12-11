## Calculate logintlik via laplace method. Bayes factor is based on
## the difference between two logintlik.

logintlik <- function(model = "linkinteraction.jag", 
                      paraest = c("beta0", "lambda1", 
                                  "lambda2", 
                                  "beta1", "beta2", "beta3", 
                                  "Inv_sig_e1", "phi"),
                      data =  c(sureriskreg, list(t=1)),
                      priorphi = "dunif"){
  m1 <- jags.model(model, data, n.chains=3)
  update(m1, 1000)
  m1res <- coda.samples(m1, paraest, n.iter=5000)
  ## posterior medians
  m1med <- summary(m1res)$quantiles[,3]
  ## posterior covariances
  m1cov <- cov(rbind(m1res[[1]], m1res[[2]], m1res[[3]]))
  ## all possible parameter names in bayeslv.Rnw. 
  paranam = c("lambda1", "lambda2",
              "int1","int2", "delta1", "delta2", "delta3",
              "beta0", "beta1", "beta2", "beta3",
              "b0", "b1", "b2", "b3",
              "g1", "g2", "g3", "g4",
              "Inv_cov1", "Inv_cov2", "phi", 
              "Inv_sig_e1", "Inv_sig_e2", "Inv_sig_ee", "Inv_sig_f1")
  pname <- names(m1med)
  rnums <- sapply(paranam, grep, pname)
  
  ## remove equal parameter Inv_cov[1:2]=Inv_cov[2:1] and
  ## fixed parameters
  if (priorphi == "dwish") {
    if (is.na(rnums$Inv_cov2[3])) {
      m1cov <- m1cov[-c(rnums$Inv_cov1[3], 
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                   m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)],
                     rnums$delta1[which(m1med[rnums$delta1]==0|
                                        m1med[rnums$delta1]==1)],
                     rnums$delta2[which(m1med[rnums$delta2]==0|
                                        m1med[rnums$delta2]==1)]),
                    -c(rnums$Inv_cov1[3],
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)],
                     rnums$delta1[which(m1med[rnums$delta1]==0|
                                        m1med[rnums$delta1]==1)],
                     rnums$delta2[which(m1med[rnums$delta2]==0|
                                        m1med[rnums$delta2]==1)])]
    } else {
     m1cov <- m1cov[-c(rnums$Inv_cov1[3], rnums$Inv_cov2[3],
                    rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                        m1med[rnums$lambda1]==1)],
                    rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                        m1med[rnums$lambda2]==1)],
                    rnums$delta1[which(m1med[rnums$delta1]==0|
                                       m1med[rnums$delta1]==1)],
                    rnums$delta2[which(m1med[rnums$delta2]==0|
                                       m1med[rnums$delta2]==1)]),
                    -c(rnums$Inv_cov1[3], rnums$Inv_cov2[3],
                    rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                        m1med[rnums$lambda1]==1)],
                    rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                        m1med[rnums$lambda2]==1)],
                    rnums$delta1[which(m1med[rnums$delta1]==0|
                                       m1med[rnums$delta1]==1)],
                    rnums$delta2[which(m1med[rnums$delta2]==0|
                                       m1med[rnums$delta2]==1)])]
    }
  } 
  if (priorphi == "gamma") {
      m1cov <- m1cov
  }
  if (priorphi == "dunif") {
     m1cov <- m1cov[-c(
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)],
                     rnums$delta1[which(m1med[rnums$delta1]==0|
                                        m1med[rnums$delta1]==1)],
                     rnums$delta2[which(m1med[rnums$delta2]==0|
                                        m1med[rnums$delta2]==1)]),
                    -c(
                     rnums$lambda1[which(m1med[rnums$lambda1]==0|
                                         m1med[rnums$lambda1]==1)],
                     rnums$lambda2[which(m1med[rnums$lambda2]==0|
                                         m1med[rnums$lambda2]==1)],
                     rnums$delta1[which(m1med[rnums$delta1]==0|
                                        m1med[rnums$delta1]==1)],
                     rnums$delta2[which(m1med[rnums$delta2]==0|
                                        m1med[rnums$delta2]==1)])]
  }

  
  ## log(f(theta*))
  logfthet <- sum(dnorm(m1med[c(
              rnums$lambda1[which(m1med[rnums$lambda1]!=0 &
                                  m1med[rnums$lambda1]!=1)],
              rnums$lambda2[which(m1med[rnums$lambda2]!=0 &
                                  m1med[rnums$lambda2]!=1)],
              rnums$delta1[which(m1med[rnums$delta1]!=0 &
                                 m1med[rnums$delta1]!=1)],
              rnums$delta2[which(m1med[rnums$delta2]!=0 &
                                 m1med[rnums$delta2]!=1)],
              rnums$int1, rnums$int2, rnums$delta3,        
              rnums$beta0, rnums$beta1, rnums$beta2,
              rnums$beta3, rnums$b0, rnums$b1, rnums$b2, rnums$b3,
              rnums$g1, rnums$g2, rnums$g3, rnums$g4)],
              0, sqrt(1/.001),
              log=TRUE)) +
              ifelse(is.na(rnums$Inv_cov1[3]),
             0, log(dwish(matrix(c(m1med[rnums$Inv_cov1]), 2, 2),2, 
                  matrix(c(1,0,0,1),2,2)))) +
              ifelse(is.na(rnums$Inv_cov2[3]),
             0, log(dwish(matrix(c(m1med[rnums$Inv_cov2]), 2, 2),2, 
                  matrix(c(1,0,0,1),2,2)))) +
          sum(dunif(m1med[rnums$phi], -1, 1, log=TRUE)) +
         sum(dgamma(m1med[rnums$Inv_sig_f1], .01, .01, log=TRUE)) + 
          sum(dgamma(m1med[rnums$Inv_sig_e1], .01, .01,log=TRUE)) +
          sum(dgamma(m1med[rnums$Inv_sig_e2], .01, .01,log=TRUE)) +
          sum(dgamma(m1med[rnums$Inv_sig_ee], .01, .01,log=TRUE)) 
  
 

  ## log-conditional likelihood, with thetas integrated out
  if (is.na(rnums$Inv_cov2[3])) {
    lcl <- 0
    if (priorphi == "dwish"){
      lambda <- matrix(c(m1med[rnums$lambda1],
                          m1med[rnums$lambda2]), 10, 2)
      phi <- solve(matrix(c(m1med[rnums$Inv_cov1]), 2, 2))
    }
    if (priorphi == "gamma"){
      lambda <- matrix(c(m1med[rnums$lambda1]), 10, 1)
      phi <- 1/m1med[rnums$Inv_sig_f1]
    }
    if (priorphi == "dunif"){
      lambda <- matrix(c(m1med[rnums$lambda1],
                          m1med[rnums$lambda2]), 10, 2)
      phi <- matrix(c(1, m1med[rnums$phi], m1med[rnums$phi], 1),
                    2, 2)
    }
    if (is.na(rnums$Inv_sig_ee[1])){
      psi <- diag(1/m1med[rnums$Inv_sig_e1])
      covmat <- lambda %*% phi %*% t(lambda) + psi
    } else {
      psi <- diag(1/m1med[rnums$Inv_sig_ee])
      covmat <- lambda %*% phi %*% t(lambda)
    }
  
    for(i in 1:data$N){
      if (!is.na(rnums$int1[1])){
        mnvec <- m1med[rnums$int1]
        lcl <- lcl + dmvnorm(data$y[i,], mnvec, covmat, log=TRUE)
      } 
      if (model == "linkinteraction.jag"){
        mnvec <- m1med[rnums$beta0] +
                 c(rep(m1med[rnums$beta1[1]],5),
                 rep(m1med[rnums$beta1[2]],5))*data$Frame[i] + 
                 m1med[rnums$beta2]*data$Num[i]  + 
                 if(is.na(m1med[rnums$beta3][1])){rep(0,10)}else{ 
                 (rep(c(1,0), each=5)*
                 m1med[rnums$beta3]*data$Frame[i]*data$Num[i])}
        lcl <- lcl + dmvnorm(data$y[i,], mnvec, covmat, log=TRUE)
      }
      if (model == "linkFramesure1.jag"){
        mnvec <- m1med[rnums$beta0] +
                 if(is.na(m1med[rnums$beta1][1])){rep(0,10)}else{
                 (c(rep(m1med[rnums$beta1[1]],
                 5), rep(0, 5))* data$Frame[i])} +
                 m1med[rnums$beta2]*data$Num[i]
        lcl <- lcl + dmvnorm(data$y[i,], mnvec, covmat, log=TRUE)
      }
      if (model == "linkFramerisky1.jag"){
        mnvec <- m1med[rnums$beta0] +
                 if(is.na(m1med[rnums$beta1][1])){rep(0,10)}else{
                 (c(rep(0, 5), rep(m1med[rnums$beta1[1]],
                 5))* data$Frame[i])} +
                 m1med[rnums$beta2]*data$Num[i]
        lcl <- lcl + dmvnorm(data$y[i,], mnvec, covmat, log=TRUE)
      }
      if (!is.na(m1med[rnums$Inv_sig_ee][1])){
        xcovmat <- m1med[rnums$g1]^2*covmat[6:10,6:10] + 
                   m1med[rnums$g2]^2*covmat[1:5,1:5] - 
                   m1med[rnums$g1]*m1med[rnums$g2]*
                   covmat[1:5,6:10] %*% 
                   covmat[6:10,1:5] +
                   (sum(m1med[rnums$g3])*(data$Num[i]-data$Numc))^2*
                   covmat[6:10,6:10] +
                   (sum(m1med[rnums$g4])*(data$Num[i]-data$Numc))^2*
                   covmat[1:5,1:5] - 
                   sum(m1med[rnums$g3])*(data$Num[i]-data$Numc)*
                   sum(m1med[rnums$g4])*
                   (data$Num[i]-data$Numc)*covmat[6:10,1:5] %*% 
                   covmat[1:5,6:10] + psi
  
        mnvec <-  m1med[rnums$b0] + rep(sum(m1med[rnums$b1]),5)*
                  data$Frame[i] + m1med[rnums$b2]*
                  (data$Num[i]-data$Numc) +
                  sum(m1med[rnums$b3])*data$Frame[i]*
                  (data$Num[i]-data$Numc) +
                  m1med[rnums$g1]*
                  m1med[rnums$beta0[6:10]] + m1med[rnums$g2]*
                  m1med[rnums$beta0[1:5]] + 
                  sum(m1med[rnums$g3])*m1med[rnums$beta0[6:10]]*
                  (data$Num[i]-data$Numc) + 
                  sum(m1med[rnums$g4])*m1med[rnums$beta0[1:5]]* 
                  (data$Num[i]-data$Numc)
        lcl <- lcl + dmvnorm(data$X[i,], mnvec, xcovmat, log=TRUE)
     }
   }
 } else {
  lcl1 <- 0
  lcl2 <- 0
  lambda1 <- matrix(c(m1med[rnums$lambda1],
                          m1med[rnums$lambda2]), 10, 2)
  lambda2 <- matrix(c(m1med[rnums$lambda1] +
             if(is.na(rnums$delta1[1])){rep(0,10)}else{
                 m1med[rnums$delta1]},
             m1med[rnums$lambda2] +
             if(is.na(rnums$delta2[1])){rep(0,10)}else{
                m1med[rnums$delta2]}),
             10, 2)
  phi1 <- solve(matrix(c(m1med[rnums$Inv_cov1]), 2, 2))
  phi2 <- solve(matrix(c(m1med[rnums$Inv_cov2]), 2, 2))
  psi1 <- diag(1/m1med[rnums$Inv_sig_e1])
  psi2 <- diag(1/m1med[rnums$Inv_sig_e2])
  covmat1 <- lambda1 %*% phi1 %*% t(lambda1) + psi1
  covmat2 <- lambda2 %*% phi2 %*% t(lambda2) + psi2
  
  for(i in 1:(data$N/2-1)){
    mnvec <-  m1med[rnums$int1]
    lcl1 <- lcl1 + dmvnorm(data$y[i,], mnvec, covmat1, log=TRUE)
  }
  for(i in (data$N/2):(data$N)){
    if(is.na(rnums$delta3[1]) & !is.na(rnums$int2[1])){
        mnvec <- m1med[rnums$int2]}
    if(!is.na(rnums$delta3[1]) & is.na(rnums$int2[1])){
         mnvec <- m1med[rnums$int1] + m1med[rnums$delta3]}
    if(is.na(rnums$delta3[1]) & is.na(rnums$int2[1])){
         mnvec <- m1med[rnums$int1]}
    lcl2 <- lcl2 + dmvnorm(data$y[i,], mnvec, covmat2, log=TRUE)
  }
  lcl <- lcl1 + lcl2
}
  ## Eq (7) of lewis + raftery
  logintlik <- (nrow(m1cov)/2)*log(2*pi) + 
                .5*log(det(m1cov)) + logfthet + lcl
  logintlik
}


#############################################################
### Calculate Bayes factor via path sampling
#bfcalc <- function(model, data, S=10, fname)
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
#
 
