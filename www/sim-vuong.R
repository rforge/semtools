##################
## SIMULATION 1 ##
##################

## Data-generating process
dgp1 <- function(parval, nobs) {
  ## Define data generating model, which should be distinct
  ## from the fitted models.

  ## weights from HS1939 data in lavaan, set 1 on variances of latent variables
  # visual <- c(.9, .498, .656)
  # textual <- c(.99, 1.102, .917)
  # speed <- c(.619, .731, .67)

  # covariances
  # v.t <- .459
  # v.s <- .471
  # t.s <- .283

  # variances : x1 to x9
  # c(.549, 1.134, .844, .371, .446, .356, .799, .488, .566)

  ## parval reflects model parameter that makes one of the fitted
  ## models preferable to the other?
  # order : c(x4 ~ x9)
  lambda <- cbind(c(.99, 1.102, .917,  parval, 0, 0), # textual
                  c(0, 0, 0,  .619, .731, .67)) # speed

  phi <- matrix(.283, 2, 2)
  diag(phi) <- 1

  # order : c(x4 ~ x9)
  theta <- diag(c(.371, .446, .356, .799, .488, .566))

  ## population covariance matrix
  popCov <- lambda %*% phi %*% t(lambda) + theta

  ## Return data matrix
  datmat <- rmvnorm(n=nobs, sigma=popCov)
  colnames(datmat) <- c("x1", "x2", "x3", "x4", "x5", "x6")

  list(d=as.data.frame(datmat), true_cov=popCov)
}

## Evaluate power on a single dgp() scenario
testpower1 <- function(parval, nobs,
                        nrep, alpha=0.05, R, verbose=TRUE) {
  checkCond <- function(code) {
    tryCatch(code,
             error = function(c) "ill-conditioned",
             warning = function(c) "ill-conditioned")
  }

  ## Function to obtain bootstrap BICs
  bootBICf <- function(D, idx) {
    boot1 <- checkCond(cfa(myModel1, data=D[idx,], std.lv=TRUE))
    boot2 <- checkCond(cfa(myModel2, data=D[idx,], std.lv=TRUE))
    if (all(!sapply(list(boot1, boot2),
                    identical, "ill-conditioned"))) {
      bic1 <- BIC(boot1)
      bic2 <- BIC(boot2)
      bic2 - bic1
    } else {
      NA
    }
  }

  ## define myfit1 and myfit2, which fit candidate models
  myModel1 <- ' F1 =~ x1 + x2 + x3 + x4
                  F2 =~ x4 + x5 + x6
                '
  myModel2 <- ' F1 =~ x1 + x2 + x3
                  F2 =~ x3 + x4 + x5 + x6
                '


  pval_pre <- mclapply(1:nrep, function(i) {
    d <- dgp1(parval, nobs)$d

    my1 <- checkCond(cfa(myModel1, data=d, std.lv=TRUE))
    my2 <- checkCond(cfa(myModel2, data=d, std.lv=TRUE))

    if(all(!sapply(list(my1, my2),
                   identical, "ill-conditioned"))) {
      ## code to compute test statistics
      res <- checkCond(vuongtest(my1, my2))
      ci <- checkCond(icci(my2, my1, 1 - 2*alpha))
      bootBIC <- boot(d, bootBICf, R=R, parallel = "multicore")
      bootBIC_ci <- quantile(bootBIC$t, probs=c(alpha, 1-alpha),
                             na.rm=TRUE, names=FALSE)

      ## code for NET: fit Model2 with the implied cov and mean from Model1
      ## Mean structure is not modeled, so it is useless and gives
      ## convergence problems.
      cov1 <- unclass(fitted(my1)$cov)
      nobs1 <- nobs(my1)
      my2_cov1 <- checkCond(cfa(myModel2,
                          sample.cov=cov1, #sample.mean=mean1,
                          sample.nobs=nobs1, std.lv=TRUE)) #, meanstructure=TRUE)

      if(all(!sapply(list(res, ci, my2_cov1), identical, "ill-conditioned"))) {
        list(p_omega = res$p_omega,
             p_LRT = res$p_LRT$A,
             p_LRT_B = res$p_LRT$B,
             BICpref = BIC(my2) > BIC(my1),
             # F value of NET
             NET = inspect(my2_cov1, "fit")[["chisq"]],
             BIC_Vuong_low = ci$BICci[1],
             BIC_Vuong_upp = ci$BICci[2],
             BIC_boot_low = bootBIC_ci[1],
             BIC_boot_upp = bootBIC_ci[2])
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, mc.cores = round(.75*detectCores()), mc.preschedule=FALSE)

  pval <- tryCatch(rbindlist(pval_pre),
                   error = function(e) {
                     print(pval_pre)
                     return(list(my1=my1, my2=my2, res=res,
                                 ci=ci, my2_cov1=my2_cov1))
                   })

  ## add a new column to pval:
  ## true BIC: subtract difference in df, see Eq (7) of Preacher + Merkle
  true_cov <- dgp1(parval, nobs)$true_cov
  dimnames(true_cov)  <- list(paste0("x", 1:6), paste0("x", 1:6))

  model1_tcov <- sem(myModel1, sample.cov=true_cov, sample.nobs=nobs,
                     sample.cov.rescale=FALSE, std.lv=TRUE)
  model2_tcov <- sem(myModel2, sample.cov=true_cov, sample.nobs=nobs,
                     sample.cov.rescale=FALSE, std.lv=TRUE)
  pval[, `:=`(tBIC = (BIC(model2_tcov) + model1_tcov@Fit@npar -
                          BIC(model1_tcov) - model2_tcov@Fit@npar))]

  ## 1,2. Save unconditional power of each test,
  ## 3. along with power of LRT conditioned on significant overlapping test
  ## NET : equivalence test with e=.001, 1 for non-equivalence, 0 for equivalence
  ## NETval : mean of F values
  rval <- pval[ , list(
      parval = parval, nobs = nobs,
      Overlap = mean(p_omega < alpha, na.rm = TRUE),
      LRT = mean(p_LRT < alpha, na.rm = TRUE),
      CondLRT = if (all(!(ol <- p_omega < alpha))) 0
                else mean(p_LRT[ol] < alpha, na.rm = TRUE),
      BIC = mean(BICpref, na.rm = TRUE),
      BICvar = var(BICpref, na.rm = TRUE),
      NET = mean(NET > .001, na.rm = TRUE),
      NETval = mean(NET, na.rm = TRUE),
      p_tBIC_Vuong = mean((BIC_Vuong_low < tBIC) & (tBIC < BIC_Vuong_upp), na.rm=TRUE),
      BIC_CI_avg_Vuong = mean(BIC_Vuong_upp - BIC_Vuong_low, na.rm=TRUE),
      BIC_CI_sd_Vuong = sqrt((sum((BIC_Vuong_upp - mean(BIC_Vuong_upp))^2) +
                                  sum((BIC_Vuong_low - mean(BIC_Vuong_low))^2))/
                                      (2*length(BIC_Vuong_upp))),
      p_tBIC_boot = mean((BIC_boot_low < tBIC) & (tBIC < BIC_boot_upp), na.rm=TRUE),
      BIC_CI_avg_boot = mean(BIC_boot_upp - BIC_boot_low, na.rm=TRUE),
      BIC_CI_sd_boot = sqrt((sum((BIC_boot_upp - mean(BIC_boot_upp))^2) +
                                 sum((BIC_boot_low - mean(BIC_boot_low))^2))/
                                     (2*length(BIC_boot_upp)))
  )]

  if(verbose) print(rval)

  return(rval)
  #return(list(pval=pval, rval=rval))
}

## Loop over scenarios
simulation1 <- function(parval=c(0, .25, .5), nobs=c(200, 500, 1000),
                        nrep=3000, alpha=0.05, R=1000, verbose=TRUE) {
  pars <- expand.grid(parval=parval, nobs=nobs)
  npars <- NROW(pars)
  pow <- vector(mode="list", length=npars)

  for(i in 1:npars) {
    if(verbose) print(pars[i,])
    pow[[i]] <- testpower1(parval=pars[i,"parval"], nobs=pars[i,"nobs"],
                            nrep=nrep, alpha=alpha, R=R, verbose=verbose)
  }

  rval <- rbindlist(pow)
  return(melt(rval, id.vars=c("parval", "nobs"),
              variable.name="test", value.name="power"))
}



##################
## SIMULATION 2 ##
##################

## Data generating function
dg2 <- function() {
  ## Data generating model, using simsem::model
  ## : Preacher & Merkle (2012) Model D

  ## BE: Regression coefficient matrix among endogenous factors
  BE <- matrix(0, 9, 9)
  BE[2:4, 1] <- NA # a  to  b,c,d
  BE[5:7, 2] <- NA # b  to  e,f,g
  BE[7:8, 3] <- NA # c  to  g,h
  BE[8:9, 4] <- NA # d  to  h,i
  BEval <- BE
  BEval[is.na(BEval)] <- 0.2
  modelD_BE <- bind(BE, BEval)

  ## BE2: Kris's actual model D
  BE2 <- matrix(0, 9, 9)
  BE2[2:4, 1] <- NA # a  to  b,c,d
  BE2[5:6, 2] <- NA # b  to  e,f
  BE2[6:7, 3] <- NA # c  to  f,g
  BE2[7:9, 4] <- NA # d  to  g,h,i
  BE2val <- BE2
  BE2val[is.na(BE2val)] <- 0.2
  modelD_BE2 <- bind(BE2, BE2val)

  ## PS: Residual variance-covariance matrix among endogenous factors
  PS <- diag(NA, 9)
  PSval <- diag(0.8, 9)
  PSval[1, 1] <- 1
  modelD_PS <- binds(PS, PSval)
  modelD_path <- model(PS=modelD_PS, BE=modelD_BE, modelType = "Path")
  modelD2_path <- model(PS=modelD_PS, BE=modelD_BE2, modelType="Path")
  ## checking model
  # dd <- generate(modelD_path, n=100000)
  # summary(analyze(modelD_path, dd)) # df=26, matched to the paper

  ## True covariance matrix for true BIC
  I_B_inv <- solve(diag(1, 9) - BEval)
  I_B2_inv <- solve(diag(1, 9) - BE2val)
  true_cov <- I_B_inv %*% PSval %*% t(I_B_inv)
  true_cov2 <- I_B2_inv %*% PSval %*% t(I_B2_inv)
  dimnames(true_cov) <- dimnames(true_cov2) <- list(paste0("y", 1:9), paste0("y", 1:9))
  # summary(sem(modelD_path@pt, sample.cov=true_cov, sample.nobs=200))


  ## objects of simsem::model have 'parameter table', @pt,
  ##  which can be used in 'model' argument of lavaan
  ##-------  model A  ----------##
  BEa <- matrix(0, 9, 9)
  BEa[3, 1] <- NA # a  to  c
  BEa[3, 2] <- NA # b  to  c
  BEa[4:9, 3] <- NA # c  to  d,e,f,g,h,i
  modelA_BE <- bind(BEa)

  PSa <- diag(NA, 9)
  PSa[1:2, 1:2] <- NA
  modelA_PS <- binds(PSa)
  modelA_path <- model(PS=modelA_PS, BE=modelA_BE, modelType = "Path")
  # analyze(modelA_path, dd) # df=27, matched to the paper
  modelA_pt <- modelA_path@pt

  ##-------  model B  ----------##
  BEb <- matrix(0, 9, 9)
  BEb[2, 1] <- NA # a  to  b
  BEb[c(3,4,6), 2] <- NA # b  to  c,d,f
  BEb[c(4:7,9), 3] <- NA # c  to  d,e,f,g,i
  BEb[7, 4] <- NA # d  to  g
  BEb[8, 5] <- NA # e  to  h
  BEb[9, 6] <- NA # f  to  i
  modelB_BE <- bind(BEb)

  PSb <- diag(NA, 9)
  PSb[7:9, 7:9] <- NA
  modelB_PS <- binds(PSb)
  modelB_path <- model(PS=modelB_PS, BE=modelB_BE, modelType = "Path")
  # analyze(modelB_path, dd) # df=21, matched to the paper
  modelB_pt <- modelB_path@pt

  ##-------  model C  ----------##
  BEc <- matrix(0, 9, 9)
  BEc[6, 3] <- NA # c  to  f
  BEc[7, 4] <- NA # d  to  g
  BEc[8, 5] <- NA # e  to  h
  BEc[9, 6:8] <- NA # i  from  g,h,i
  modelC_BE <- bind(BEc)

  ## Kris's real Model C
  BEc2 <- matrix(0, 9, 9)
  BEc2[6, 1] <- NA # a  to  f
  BEc2[7, 2] <- NA # b  to  g
  BEc2[8, 3] <- NA # c  to  h
  BEc2[9, 6:8] <- NA # i  from  g,h,i
  modelC2_BE <- bind(BEc2)

  PSc <- diag(NA, 9)
  PSc[1:5, 1:5] <- NA
  modelC_PS <- binds(PSc)
  modelC_path <- model(PS=modelC_PS, BE=modelC_BE, modelType = "Path")
  # analyze(modelC_path, dd) # df=20, matched to the paper
  modelC2_path <- model(PS=modelC_PS, BE=modelC2_BE, modelType="Path")
  modelC_pt <- modelC_path@pt
  modelC2_pt <- modelC2_path@pt

  list(modelA_pt = modelA_pt,
       modelB_pt = modelB_pt,
       modelC_pt = modelC_pt,
       modelC2_pt = modelC2_pt,
       modelD_path = modelD_path,
       modelD2_path = modelD2_path,
       true_cov = true_cov,
       true_cov2 = true_cov2)
}


## Evaluate power
testpower2 <- function(nrep=1000, alpha=.05, nobs=200, verbose=TRUE, R=1000, ...)
{
  dgl <- dg2()
  modelA_pt <- dgl$modelA_pt
  modelB_pt <- dgl$modelB_pt
  modelC_pt <- dgl$modelC_pt
  modelD_path <- dgl$modelD_path
  true_cov <- dgl$true_cov

  ## Function to obtain bootstrap BICs
  bootBICf <- function(D, idx) {
    a <- BIC(sem(modelA_pt, data=D[idx,]))
    b <- BIC(sem(modelB_pt, data=D[idx,]))
    c <- BIC(sem(modelC_pt, data=D[idx,]))
    ab <- b - a
    bc <- c - b
    ca <- a - c
    c(ab=ab, bc=bc, ca=ca)
  }

  ## a function checks error and warning
  check_condition <- function(code) {
    tryCatch(code,
             error = function(c) "ill-conditioned",
             warning = function(c) "ill-conditioned")
  }

  pval_pre <- mclapply(1:nrep, function(i){

    ## generating data from modelD
    d <- generate(modelD_path, n=nobs)

    ## fit candidate models
    modelA <- check_condition(sem(modelA_pt, data=d))
    modelB <- check_condition(sem(modelB_pt, data=d))
    modelC <- check_condition(sem(modelC_pt, data=d))

    if(all(!sapply(list(modelA, modelB, modelC),
                   identical, "ill-conditioned"))) {

      ## compute test statistics
      resAB <- try(vuongtest(modelB, modelA))
      resBC <- try(vuongtest(modelC, modelB))
      resCA <- try(vuongtest(modelA, modelC))

      ciAB <- try(icci(modelB, modelA, 1 - 2*alpha))
      ciBC <- try(icci(modelC, modelB, 1 - 2*alpha))
      ciCA <- try(icci(modelA, modelC, 1 - 2*alpha))

      ## code for NET: fit Model2 with the implied cov and mean from Model1
      ## Mean structure is not modeled, so it is useless and gives
      ## convergence problems.
      covA <- unclass(fitted(modelA)$cov)
      covB <- unclass(fitted(modelB)$cov)
      modelB_covA <- sem(modelB_pt, sample.cov=covA, sample.nobs=nobs)
      modelC_covA <- sem(modelC_pt, sample.cov=covA, sample.nobs=nobs)
      modelC_covB <- sem(modelC_pt, sample.cov=covB, sample.nobs=nobs)

      bootBIC <- boot(d, bootBICf, R=R, parallel = "multicore")

      ## CI, median and mean for each replication
      ## 100*(1-alpha)% two-sided CI to match one-sided test with alpha
      bootBIC_ci <- apply(bootBIC$t, 2, quantile,
                          probs=c(alpha, 0.5, 1-alpha),
                          na.rm=TRUE, names=FALSE)
      bootBIC_mean <- apply(bootBIC$t, 2, mean, na.rm=TRUE)

      ## Logic: if nothing is a try-error, retain the results.
      ## Otherwise, discard results
      if(all(!sapply(list(resAB, resBC, resCA, ciAB, ciBC, ciCA),
                     inherits, "try-error"))) {
        list(p_Omega_AB = resAB$p_omega,
             p_Omega_BC = resBC$p_omega,
             p_Omega_CA = resCA$p_omega,
             p_LRT_AB = resAB$p_LRT$B,
             p_LRT_BC = resBC$p_LRT$B,
             p_LRT_CA = resCA$p_LRT$B,
             BIC_AB = as.numeric(BIC(modelB) > BIC(modelA)),
             BIC_BC = as.numeric(BIC(modelC) > BIC(modelB)),
             BIC_CA = as.numeric(BIC(modelA) > BIC(modelC)),
             NET_BA = inspect(modelB_covA, "fit")[["chisq"]],
             NET_CA = inspect(modelC_covA, "fit")[["chisq"]],
             NET_CB = inspect(modelC_covB, "fit")[["chisq"]],
             BIC_AB_Vuong_low = ciAB$BICci[1],
             BIC_AB_Vuong_upp = ciAB$BICci[2],
             BIC_BC_Vuong_low = ciBC$BICci[1],
             BIC_BC_Vuong_upp = ciBC$BICci[2],
             BIC_CA_Vuong_low = ciCA$BICci[1],
             BIC_CA_Vuong_upp = ciCA$BICci[2],
             BIC_AB_boot_low = bootBIC_ci[1, 1],
             BIC_AB_boot_med = bootBIC_ci[2, 1],
             BIC_AB_boot_upp = bootBIC_ci[3, 1],
             BIC_AB_boot_mean = bootBIC_mean[[1]],
             BIC_BC_boot_low = bootBIC_ci[1, 2],
             BIC_BC_boot_med = bootBIC_ci[2, 2],
             BIC_BC_boot_upp = bootBIC_ci[3, 2],
             BIC_BC_boot_mean = bootBIC_mean[[2]],
             BIC_CA_boot_low = bootBIC_ci[1, 3],
             BIC_CA_boot_med = bootBIC_ci[2, 3],
             BIC_CA_boot_upp = bootBIC_ci[3, 3],
             BIC_CA_boot_mean = bootBIC_mean[[3]])
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, mc.cores = round(.75*detectCores()), mc.preschedule=FALSE)

  ## This removes iterations that resulted in a try-error:
  pval <- rbindlist(pval_pre)

  ## true BIC difference
  modelA_tcov <- sem(modelA_pt, sample.cov=true_cov, sample.nobs=nobs, sample.cov.rescale=FALSE)
  modelB_tcov <- sem(modelB_pt, sample.cov=true_cov, sample.nobs=nobs, sample.cov.rescale=FALSE)
  modelC_tcov <- sem(modelC_pt, sample.cov=true_cov, sample.nobs=nobs, sample.cov.rescale=FALSE)
  ## add three new columns to pval:
  ## true BIC: subtract difference in df, see Eq (7) of Preacher + Merkle
  ## May not be needed for vuong intervals
  pval[, `:=`(tBIC_AB=(BIC(modelB_tcov) + modelA_tcov@Fit@npar -
                       BIC(modelA_tcov) - modelB_tcov@Fit@npar),
              tBIC_BC=(BIC(modelC_tcov) + modelB_tcov@Fit@npar -
                       BIC(modelB_tcov) - modelC_tcov@Fit@npar),
              tBIC_CA=(BIC(modelA_tcov) + modelC_tcov@Fit@npar -
                       BIC(modelC_tcov) - modelA_tcov@Fit@npar))]

  ## Now summarize results across reps
  ## NB: && does not work the same as & here!!
  rval <- pval[, list(
      overlap_AB = mean(p_Omega_AB < alpha, na.rm=TRUE),
      overlap_BC = mean(p_Omega_BC < alpha, na.rm=TRUE),
      overlap_CA = mean(p_Omega_CA < alpha, na.rm=TRUE),
      lrt_AB = mean(p_LRT_AB < alpha, na.rm=TRUE),
      lrt_BC = mean(p_LRT_BC < alpha, na.rm=TRUE),
      lrt_CA = mean(p_LRT_CA < alpha, na.rm=TRUE),
      BIC_AB = mean(BIC_AB, na.rm=TRUE),
      BIC_BC = mean(BIC_BC, na.rm=TRUE),
      BIC_CA = mean(BIC_CA, na.rm=TRUE),
      NET_BA = mean(NET_BA > .001, na.rm=TRUE),
      NET_CA = mean(NET_CA > .001, na.rm=TRUE),
      NET_CB = mean(NET_CB > .001, na.rm=TRUE),
      p_tBIC_AB_Vuong = mean((BIC_AB_Vuong_low < tBIC_AB) & (tBIC_AB < BIC_AB_Vuong_upp), na.rm=TRUE),
      p_tBIC_BC_Vuong = mean((BIC_BC_Vuong_low < tBIC_BC) & (tBIC_BC < BIC_BC_Vuong_upp), na.rm=TRUE),
      p_tBIC_CA_Vuong = mean((BIC_CA_Vuong_low < tBIC_CA) & (tBIC_CA < BIC_CA_Vuong_upp), na.rm=TRUE),
      BIC_AB_CI_avg_Vuong = mean(BIC_AB_Vuong_upp - BIC_AB_Vuong_low, na.rm=TRUE),
      BIC_BC_CI_avg_Vuong = mean(BIC_BC_Vuong_upp - BIC_BC_Vuong_low, na.rm=TRUE),
      BIC_CA_CI_avg_Vuong = mean(BIC_CA_Vuong_upp - BIC_CA_Vuong_low, na.rm=TRUE),
      BIC_AB_CI_sd_Vuong = sqrt((sum((BIC_AB_Vuong_upp - mean(BIC_AB_Vuong_upp))^2) +
                                     sum((BIC_AB_Vuong_low - mean(BIC_AB_Vuong_low))^2))/
                                         (2*length(BIC_AB_Vuong_upp))),
      BIC_BC_CI_sd_Vuong = sqrt((sum((BIC_BC_Vuong_upp - mean(BIC_BC_Vuong_upp))^2) +
                                     sum((BIC_BC_Vuong_low - mean(BIC_BC_Vuong_low))^2))/
                                         (2*length(BIC_BC_Vuong_upp))),
      BIC_CA_CI_sd_Vuong = sqrt((sum((BIC_CA_Vuong_upp - mean(BIC_CA_Vuong_upp))^2) +
                                     sum((BIC_CA_Vuong_low - mean(BIC_CA_Vuong_low))^2))/
                                         (2*length(BIC_CA_Vuong_upp))),
      p_tBIC_AB_boot = mean((BIC_AB_boot_low < tBIC_AB) & (tBIC_AB < BIC_AB_boot_upp), na.rm=TRUE),
      p_tBIC_BC_boot = mean((BIC_BC_boot_low < tBIC_BC) & (tBIC_BC < BIC_BC_boot_upp), na.rm=TRUE),
      p_tBIC_CA_boot = mean((BIC_CA_boot_low < tBIC_CA) & (tBIC_CA < BIC_CA_boot_upp), na.rm=TRUE),
      BIC_AB_CI_avg_boot = mean(BIC_AB_boot_upp - BIC_AB_boot_low, na.rm=TRUE),
      BIC_BC_CI_avg_boot = mean(BIC_BC_boot_upp - BIC_BC_boot_low, na.rm=TRUE),
      BIC_CA_CI_avg_boot = mean(BIC_CA_boot_upp - BIC_CA_boot_low, na.rm=TRUE),
      BIC_AB_CI_sd_boot = sqrt((sum((BIC_AB_boot_upp - mean(BIC_AB_boot_upp))^2) +
                                    sum((BIC_AB_boot_low - mean(BIC_AB_boot_low))^2))/
                                        (2*length(BIC_AB_boot_upp))),
      BIC_BC_CI_sd_boot = sqrt((sum((BIC_BC_boot_upp - mean(BIC_BC_boot_upp))^2) +
                                    sum((BIC_BC_boot_low - mean(BIC_BC_boot_low))^2))/
                                        (2*length(BIC_BC_boot_upp))),
      BIC_CA_CI_sd_boot = sqrt((sum((BIC_CA_boot_upp - mean(BIC_CA_boot_upp))^2) +
                                    sum((BIC_CA_boot_low - mean(BIC_CA_boot_low))^2))/
                                        (2*length(BIC_CA_boot_upp)))
  )]

  if(verbose) print(rval)

  return(as.numeric(rval))

  ## return(list(rval = rval,
  ##             boot = pval[, list(AB_median = BIC_AB_boot_med,
  ##                                AB_mean = BIC_AB_boot_mean,
  ##                                BC_median = BIC_BC_boot_med,
  ##                                BC_mean = BIC_BC_boot_mean,
  ##                                CA_median = BIC_CA_boot_med,
  ##                                CA_mean = BIC_CA_boot_mean,
  ##                                N = nobs)]) )
}



## -----------------------------------------------
## Loop over scenarios
simulation2 <- function(nobs = c(200, 500, 1000), nrep=1000, verbose = TRUE, R=1000, ...)
{
  prs <- expand.grid(nobs = nobs)
  nprs <- nrow(prs)

  ## Define which test results to save, plus entries of pow[[i]] reflecting those
  ## tests
  test <- c("overlap_AB", "overlap_BC", "overlap_CA", "lrt_AB", "lrt_BC",
            "lrt_CA", "NET_BA", "NET_CA", "NET_CB",
            "p_tBIC_AB_Vuong", "p_tBIC_BC_Vuong", "p_tBIC_CA_Vuong", "BIC_AB_CI_avg_Vuong",
            "BIC_BC_CI_avg_Vuong", "BIC_CA_CI_avg_Vuong", "BIC_AB_CI_sd_Vuong",
            "BIC_BC_CI_sd_Vuong", "BIC_CA_CI_sd_Vuong", "p_tBIC_AB_boot",
            "p_tBIC_BC_boot", "p_tBIC_CA_boot", "BIC_AB_CI_avg_boot", "BIC_BC_CI_avg_boot",
            "BIC_CA_CI_avg_boot", "BIC_AB_CI_sd_boot", "BIC_BC_CI_sd_boot",
            "BIC_CA_CI_sd_boot")
  wentries <- c(1:6, 10:30)
  ntest <- length(test)

  ## Parallelization now occurs inside testpower2
  pow <- vector(mode="list", length=nprs)
  for(i in 1:nprs) {
    if(verbose) print(prs[i,])
    pow[[i]] <- testpower2(nobs = prs$nobs[i], nrep=nrep, verbose = verbose, R=R, ...)
  }

  rval <- data.frame(nobs = rep(as.numeric(prs$nobs), each = ntest))
  rval$test <- factor(rep(test, nprs), levels = test)
  rval$nobs <- factor(rval$nobs)
  rval$result <- unlist(lapply(pow, function(x) x[wentries]))
  return(rval)
  ##return(rbindlist(pow))

  ## return(list(pow = rbindlist(lapply(pow, `[[`, "rval")),
  ##             boot = rbindlist(lapply(pow, `[[`, "boot"))) )

}



##################
## SIMULATION 3 ##
##################

## Data generating function
dg3 <- function(parval) {
  ## Data generating model, using simsem::model
  ## : Preacher & Merkle (2012) Model B
  ## : varing covariances between y7, y8, y9

  ## BE: Regression coefficient matrix among endogenous factors
  BE <- matrix(0, 9, 9)
  BE[2, 1] <- NA # a  to  b
  BE[c(3,4,6), 2] <- NA # b  to  c,d,f
  BE[c(4:7,9), 3] <- NA # c  to  d,e,f,g,i
  BE[7, 4] <- NA # d  to  g
  BE[8, 5] <- NA # e  to  h
  BE[9, 6] <- NA # f  to  i
  BEval <- BE
  BEval[is.na(BEval)] <- 0.2
  model_BE <- bind(BE, BEval)

  ## PS: Residual variance-covariance matrix among endogenous factors
  PS <- diag(NA, 9)
  PS[7:9, 7:9] <- NA
  PSval <- diag(0.8, 9)
  PSval[matrix(c(7,8,7,9, 8,7,8,9, 9,7,9,8), ncol=2, byrow=TRUE)] <- parval
  model_PS <- binds(PS, PSval)
  model_path <- model(PS=model_PS, BE=model_BE, modelType = "Path")

  ## True covariance matrix for true BIC
  I_B_inv <- chol2inv(chol(diag(1, 9) - BEval))
  true_cov <- I_B_inv %*% PSval %*% t(I_B_inv)
  dimnames(true_cov) <- list(paste0("y", 1:9), paste0("y", 1:9))

  ## objects of simsem::model have 'parameter table', @pt,
  ##  which can be used in 'model' argument of lavaan
  ##-----  model B  -----##
  modelB_BE <- bind(BE)
  modelB_PS <- binds(PS)
  modelB_path <- model(PS=modelB_PS, BE=modelB_BE, modelType = "Path")
  modelB_pt <- modelB_path@pt

  ##---  model B without covariances between y7 ~ y9 (nested)  ----##
  modelBn_BE <- modelB_BE
  PSbn <- diag(NA, 9)
  modelBn_PS <- binds(PSbn)
  modelBn_path <- model(PS=modelBn_PS, BE=modelBn_BE, modelType = "Path")
  modelBn_pt <- modelBn_path@pt

  list(modelB_pt = modelB_pt,
       modelBn_pt = modelBn_pt,
       model_path = model_path,
       true_cov = true_cov)
}


## Evaluate power
testpower3 <- function(parval=0, nobs=200,
                       nrep=1000, alpha=.05, verbose=TRUE, mD=FALSE, ...)
{
  dgl <- dg3(parval)
  modelB_pt <- dgl$modelB_pt
  modelBn_pt <- dgl$modelBn_pt
  model_path <- dgl$model_path
  true_cov <- dgl$true_cov
  if(mD){
    ## If mD, generate data from Model D, sim2.
    ## Now the full model is incorrect.
    dgld <- dg2()
    model_path <- dgld$modelD_path
    true_cov <- dgld$true_cov
  }

  ## a function checks error and warning
  check_condition <- function(code) {
    tryCatch(code,
             error = function(c) "ill-conditioned",
             warning = function(c) "ill-conditioned")
  }

  pval_pre <- mclapply(1:nrep, function(i) {
    ## generating data from modelD
    d <- generate(model_path, n=nobs)

    ## fit candidate models
    modelB <- check_condition(sem(modelB_pt, data=d))
    modelBn <- check_condition(sem(modelBn_pt, data=d))

    if(all(!sapply(list(modelB, modelBn),
                   identical, "ill-conditioned"))) {
      ## compute test statistics
      res <- try(vuongtest(modelBn, modelB, nested=TRUE))
      p_Chidiff <- pchisq(fitMeasures(modelBn)[["chisq"]] -
                              fitMeasures(modelB)[["chisq"]],
                          df=3, lower.tail=FALSE)

      ## Logic: if nothing is a try-error, retain the results.
      ## Otherwise, discard results
      if(all(!sapply(list(res), inherits, "try-error"))) {
        list(p_Omega = res$p_omega,
             p_LRT = res$p_LRT$A,
             p_Chidiff = p_Chidiff,
             BIC = as.numeric(BIC(modelBn) > BIC(modelB)))
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, mc.cores = round(.75*detectCores()), mc.preschedule=FALSE)

  ## This removes iterations that resulted in a try-error:
  pval <- rbindlist(pval_pre)

  ## true BIC difference
  modelB_tcov <- sem(modelB_pt, sample.cov=true_cov, sample.nobs=nobs, sample.cov.rescale=FALSE)
  modelBn_tcov <- sem(modelBn_pt, sample.cov=true_cov, sample.nobs=nobs, sample.cov.rescale=FALSE)

  ## Now summarize results across reps
  ## NB: && does not work the same as & here!!
  rval <- pval[, list(parval = parval, nobs = nobs,
                      Overlap = mean(p_Omega < alpha, na.rm=TRUE),
                      Chidiff = mean(p_Chidiff < alpha, na.rm=TRUE),
                      LRT = mean(p_LRT < alpha, na.rm=TRUE),
                      BIC = mean(BIC, na.rm=TRUE)
                      )]

  if(verbose) print(rval)

  return(rval)
}



## -----------------------------------------------
## Loop over scenarios
simulation3 <- function(parval=0, nobs=200,
                        nrep=1000, alpha=0.05, verbose = TRUE, mD = FALSE, ...)
{
  pars <- expand.grid(parval=parval, nobs=nobs)
  npars <- NROW(pars)
  pow <- vector(mode="list", length=npars)

  for(i in 1:npars) {
    #if(verbose) print(pars[i,])
    pow[[i]] <- testpower3(parval=pars[i,"parval"], nobs=pars[i,"nobs"],
                           nrep=nrep, alpha=alpha, verbose=verbose, mD = mD, ...)
  }

  rval <- rbindlist(pow)
  return(melt(rval, id.vars=c("parval", "nobs"),
              variable.name="test", value.name="power"))
}
