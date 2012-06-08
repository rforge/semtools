
## Data-generating process
dgp <- function(nobs = 200, diff = 3, nlevels=10)
{
  # Generates data from a factor analysis model that violates
  # measurement invariance.  Also generates an ordinal auxiliary
  # variable related to the invariance (loosely called "age").
  stopifnot(require("mvtnorm"))
  
  # ses is SE increase in parameters for each increasing value
  # of "age" past the halfway point.
  sampsize <- nobs
  ses <- diff
  sampsize.per.nobs <- nobs %/% nlevels

  ## Deal with cases where nobs/nlevels is not an integer
  age <- rep(NA, sampsize)
  age[1:(nlevels*sampsize.per.nobs)] <- rep(1:nlevels, each=sampsize.per.nobs)
  age[is.na(age)] <- sample((1:nlevels),sum(is.na(age)))
  half.level <- ifelse(nlevels %% 2 == 1, (nlevels %/% 2)+1, nlevels/2)
  
  # Define parameter vectors/matrices:
  mu <- c(29.32,24.70,14.84,10.59,19.30,18.01)

  # Loadings:
  lambda <- matrix(0,6,2)
  lambda[,1] <- c(4.92,2.96,5.96,0,0,0)
  lambda[,2] <- c(0,0,0,3.24,4.32,7.21)

  phi <- matrix(.48,2,2)
  diag(phi) <- 1

  ## Base variances
  psi <- matrix(0,6,6)
  diag(psi) <- c(26.77,13.01,30.93,3.17,8.82,22.5)

  ## Initialize data matrix:
  datmat <- matrix(NA,sampsize,6)
  colnames(datmat) <- c("x1", "x2", "x3", "y1", "y2", "y3")

  ## Generate z and u vectors for invariant individuals:
  n.invariant <- sum(age <= half.level)
  z <- t(rmvnorm(n.invariant,rep(0,2),phi))
  u <- t(rmvnorm(n.invariant,rep(0,6),psi))

  ## Based on z and u, generate data:
  tmp.ind <- which(age <= half.level)
  for (i in 1:n.invariant){
    datmat[tmp.ind[i],] <- mu + lambda%*%z[,i] + u[,i]
  }

  ## MODEL CODE FOR FINDING APPROPRIATE SES
  ##lambda.o[,1] <- lambda.y[,1] - (ses * c(8.66, 5.52, 9.19, rep(0, 3))/sqrt(nobs/2))
  ## Now do the same for invariant-violating levels
  n.vary <- nlevels - (half.level + 1)
  each.vary <- ses/n.vary
  for (i in (half.level+1):nlevels){
    tmp.n <- sum(age==i)
    tmp.ind <- which(age==i)
    z <- t(rmvnorm(tmp.n,rep(0,2),phi))

    ## SES found by fitting model to population S, multiplying
    ## observed SES by sqrt(n)
    ## FIXME: Should ses refer to maximum difference in 'end' category, or to se increase
    ##        with each category (currently coded as difference in end category)?
    diag(psi) <- diag(psi) + (each.vary*c(62.1,26.6,82.1,9.05,19.35,51.59))/sqrt(tmp.n)
    u <- t(rmvnorm(tmp.n,rep(0,6),psi))
    for (j in 1:tmp.n){
      datmat[tmp.ind[j],] <- mu + lambda%*%z[,j] + u[,j]
    }
  }
  
  cbind(as.data.frame(datmat), age)
}

## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 5000, size = 0.05, ordfun = NULL, verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, 3 * nrep), ncol = 3)
  colnames(pval) <- c("ordmax","catdiff","lrt")
  
  for(i in 1:nrep) {
    d <- dgp(...)

    mz <- ordfit(d)

    ## Requires estfun.lavaan
    ord_gefp  <- try(gefp(mz, fit = NULL, vcov = info.mzfit, order.by = d$age, sandwich = FALSE, parm = 7:12), silent = TRUE)

    if(!inherits(ord_gefp, "try-error")) {
      pval[i, 1] <- sctest(ord_gefp,  functional = ordfun)$p.value
      pval[i, 2] <- sctest(ord_gefp,  functional = catL2BB(ord_gefp))$p.value
      pval[i, 3] <- mz$lrt.p
    }
  }
  rval <- colMeans(pval < size, na.rm = TRUE)
  if(verbose) print(rval)
  
  return(rval)
}

## Loop over scenarios
simulation <- function(diff = seq(0, 5, by = 0.5),
  nobs = c(100, 500, 1000), nlevels = c(5, 10, 15), verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs)
  nprs <- nrow(prs)
  
  test <- c("ordmax","catdiff","lrt")
  ntest <- length(test)

  pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
  tmp.level <- 0
  for(i in 1:nprs) {
    if(verbose) print(prs[i,])

    ## Only simulate new critical values if you need it.
    if (tmp.level != prs$nlevels[i]){
      tmpdat <- dgp(prs$nobs[i], prs$diff[i], prs$nlevels[i])
      tmpmod <- ordfit(tmpdat)
      tmpgefp <- gefp(tmpmod, fit = NULL, vcov = info.mzfit, order.by = tmpdat$age, sandwich = FALSE, parm = 7:12)

      critvals <- ordL2BB(tmpgefp, nobs = 1000)
    }
    tmp.level <- prs$nlevels[i]

    pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i], nlevels = prs$nlevels[i], ordfun = critvals, verbose = verbose, ...)
  }

  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(test, each = nprs))
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
source("../www/estfun.lavaan.R")
source("../prelim_code/efpFunctional-cat.R")

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
