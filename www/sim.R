## Data-generating process
dgp <- function(nobs = 200, diff = 3)
{
  # Generates data from a factor analysis model that violates
  # measurement invariance.  Also generates a continuous auxiliary
  # variable related to the invariance (loosely called "age").
  
  # ses is number of SEs between group 1's parameters and group 2's
  # parameters.
  sampsize <- nobs
  ses <- diff
  
  # Could eventually manipulate this so break point is not directly
  # in the middle
  num.y <- round(sampsize/2, 0)
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
  #old# lambda.o[,1] <- lambda.y[,1] - (ses * c(.87,.6,1.32,rep(0,3)))
  lambda.o[,1] <- lambda.y[,1] - (ses * c(8.66, 5.52, 9.19, rep(0, 3))/sqrt(nobs/2))
  lambda.o[,2] <- c(0,0,0,3.24,4.32,7.21)
  
  phi <- matrix(.48,2,2)
  diag(phi) <- 1

  psi <- matrix(0,6,6)
  diag(psi) <- c(26.77,13.01,30.93,3.17,8.82,22.5)

  # Initialize data matrix:
  datmat <- matrix(0,sampsize,6)
  colnames(datmat) <- c("x1", "x2", "x3", "y1", "y2", "y3")

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
  
  cbind(as.data.frame(datmat), age)
}

## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 5000, size = 0.05, verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, 9 * nrep), ncol = 9)
  colnames(pval) <- c("dmax3", "CvM3", "supLM3", "dmax13", "CvM13", "supLM13", "dmax19", "CvM19", "supLM19")
  
  for(i in 1:nrep) {
    d <- dgp(...)
    mz <- mzfit(d)
    mz_gefp3  <- try(gefp(mz, fit = NULL, vcov = info.mzfit, order.by = d$age, sandwich = FALSE, parm = 1:3), silent = TRUE)
    mz_gefp13 <- try(gefp(mz, fit = NULL, vcov = info.mzfit, order.by = d$age, sandwich = FALSE, parm = 1:13), silent = TRUE)
    mz_gefp19 <- try(gefp(mz, fit = NULL, vcov = info.mzfit, order.by = d$age, sandwich = FALSE), silent = TRUE)

    if(!inherits(mz_gefp3, "try-error")) {
      pval[i, 1] <- sctest(mz_gefp3,  functional = maxBB)$p.value
      pval[i, 2] <- sctest(mz_gefp3,  functional = meanL2BB)$p.value
      pval[i, 3] <- sctest(mz_gefp3,  functional = supLM(0.1))$p.value
    }
    if(!inherits(mz_gefp13, "try-error")) {
      pval[i, 4] <- sctest(mz_gefp13,  functional = maxBB)$p.value
      pval[i, 5] <- sctest(mz_gefp13,  functional = meanL2BB)$p.value
      pval[i, 6] <- sctest(mz_gefp13,  functional = supLM(0.1))$p.value
    }
    if(!inherits(mz_gefp19, "try-error")) {
      pval[i, 7] <- sctest(mz_gefp19, functional = maxBB)$p.value
      pval[i, 8] <- sctest(mz_gefp19, functional = meanL2BB)$p.value
      pval[i, 9] <- sctest(mz_gefp19, functional = supLM(0.1))$p.value
    }
  }
  rval <- colMeans(pval < size, na.rm = TRUE)
  if(verbose) print(rval)
  
  return(rval)
}

## Loop over scenarios
simulation <- function(diff = seq(0, 4, by = 0.25),
  nobs = c(50, 100, 200, 500), verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nobs = nobs)
  nprs <- nrow(prs)
  
  test <- c("dmax3", "CvM3", "supLM3", "dmax13", "CvM13", "supLM13", "dmax19", "CvM19", "supLM19")
  ntest <- length(test)

  pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  for(i in 1:nprs) {
    if(verbose) print(prs[i,])
    pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i], verbose = verbose, ...)
  }

  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(c("dmax", "CvM", "supLM", "dmax", "CvM", "supLM", "dmax", "CvM", "supLM"), each = nprs),
    levels = c("dmax", "CvM", "supLM"))
  rval$pars <- factor(rep(c("3", "3", "3", "13", "13", "13", "19", "19", "19"), each = nprs), levels = c("3", "13", "19"))
  rval$nobs <- factor(rval$nobs)
  rval$power <- as.vector(pow)
  return(rval)
}

if(FALSE) {
## code and packages
library("lavaan")
library("strucchange")
library("mvtnorm")
source("mz.R")

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
mz_sim <- simulation()
save(mz_sim, file = "mz_sim.rda")
}

if(FALSE) {
## load simulation results
load("mz_sim.rda")

## display result in table
## all without rounding
ftable(xtabs(power ~ nobs + pars + test + diff, data = mz_sim), col.vars = "diff")
## subset with rounding, in percent
round(ftable(100 * xtabs(power ~ nobs + pars + test + diff,
  data = mz_sim, subset = diff %in% c(seq(0, 3.5, by = 0.5))),
  col.vars = "diff"), digits = 1)

## display result in graphic
library("lattice")
trellis.par.set(theme = canonical.theme(color = FALSE))
## all
xyplot(power ~ diff | pars + nobs, group = ~ test, data = mz_sim, type = "b",
       xlab="Violation Magnitude", ylab="Power")
}
