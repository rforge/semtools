## Data-generating process
dgp1 <- function(nobs = 9600, diff = 3, nlevels = 8, parms="gamma", meanch=0)
{
  ## Generates data from a factor analysis model that violates
  ## measurement invariance.  Also generates an ordinal auxiliary
  ## variable related to the invariance (loosely called "age").

  # require this package to use inv.logit function
  stopifnot(require("boot"))

  ses <- diff
  sampsize <- nobs
  sampsize.per.nobs <- nobs %/% nlevels

  ## Deal with cases where nobs/nlevels is not an integer
  age <- rep(NA, sampsize)
  age[1:(nlevels*sampsize.per.nobs)] <- rep(1:nlevels,
                                            each=sampsize.per.nobs)
  if (any(is.na(age))) age[is.na(age)] <- sample((1:nlevels),
                                                   sum(is.na(age)))
  half.level <- ifelse(nlevels %% 2 == 1, (nlevels %/% 2)+1,
                                                         nlevels/2)

  ## Define parameter vectors/matrices:
  theta <- c(2, 1.5, 1, 0.5, 0.3,
             -2, -1, 0, 1, 2)
  if (parms=="alphaall"|parms=="gammaall"){
  # equal difference across items. This difference number is tricky.
      asymp <- rep(2.459675, 10)
      } else {
      asymp <- c(13.640015, 7.602631, 4.472136, 3.577709,
                 4.472136, 8.944272, 3.801316, 2.459675,
                 2.459675, 3.130495)
  }


  # alpha:
  alpha <- theta[1:5]

  # gamma
  gamma <- theta[6:10]


  ## Initialize data matrix:
  datmat <- matrix(NA,sampsize,5)
  colnames(datmat) <- c("x1", "x2", "x3", "x4", "x5")

  ## Generate data for invariant individuals:
  n.invariant <- sum(age <= half.level)
  eta <- rnorm(n.invariant, 0, 1)
  prob.ind <- matrix(NA, nrow = n.invariant, ncol = 5)

  ## Based on eta, alpha, gamma, generate data:
  tmp.ind <- which(age <= half.level)
  for (i in 1:n.invariant){
    prob.ind[i,] <- inv.logit(gamma + eta[i]*alpha)
    for (j in 1:5){
      datmat[tmp.ind[i],j] <- rbinom(1, 1, prob.ind[i,j])
    }
  }

  ## Now do the same for invariant-violating levels
  n.vary <- nlevels - half.level

  for (i in (half.level+1):nlevels){
    tmp.n <- sum(age==i)
    tmp.ind <- which(age==i)


    if (parms=="alpha"){parmsvector=c(rep(0,2), 1, rep(0,7))}
    if (parms=="gamma"){parmsvector=c(rep(0,7), 1, rep(0,2))}
    if (parms=="gammaall"){parmsvector=c(rep(0,5), rep(1,5))}
    if (parms=="alphaall"){parmsvector=c(rep(1,5), rep(0,5))}
    if (parms=="mean"){parmsvector=rep(0,10)}

    theta.tmp <- theta + (ses*asymp*parmsvector/sqrt(tmp.n))
    # Differentiate and difficulty:
    tmp.alpha <- theta.tmp[1:5]
    tmp.gamma <- theta.tmp[6:10]

    if (meanch==3){eta.tmp <- rnorm(tmp.n, -ses*0.5, 1)}
    if (meanch==2){eta.tmp <- rnorm(tmp.n, -1, 2)}
    if (meanch==1){eta.tmp <- rnorm(tmp.n, -1, 1)}
    if (meanch==0){eta.tmp <- rnorm(tmp.n, 0, 1)}
    prob.tmp <- matrix(NA, nrow = tmp.n, ncol = 5)

    for (j in 1:tmp.n){
       prob.tmp[j,] <- inv.logit(tmp.gamma + eta.tmp[j]*tmp.alpha)
      for (k in 1:5){
      datmat[tmp.ind[j],k] <- rbinom(1, 1, prob.tmp[j, k])
    }
  }
}
  cbind(as.data.frame(datmat), age)
}



## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 20, size = 0.05, ordfun = NULL, parnum = parnum,
                      test = test, estimator = estimator, meantest = meantest, meanch = meanch,
                      sim ="sim1", verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, length(test) * nrep), ncol = length(test))
  colnames(pval) <- test
  dfun <- paste("dgp", substr(sim,4,4), sep="")

  for(i in 1:nrep) {
    d <- do.call(dfun, list(...))
    if (estimator == "MML"){
      mz <- tryCatch(try(mzirtfit.ltm(d)), error=function(e) e, warning=function(w) w)
    }
    if (estimator == "PML" & meantest == 0){
      mz <- tryCatch(try(mzirtfit.pml(d)), error=function(e) e, warning=function(w) w)
    }
    if (estimator == "PML" & meantest == 1){
      mz <- tryCatch(try(mzgroupmeanirtfit.pml(d)), error=function(e) e, warning=function(w) w)
    }
    if (estimator == "PML" & meantest == 2){
      mz <- tryCatch(try(mzgroupirtfit.pml(d)), error=function(e) e, warning=function(w) w)
    }
    if (estimator == "PML" & meantest == 3){
      mz <-tryCatch(try(mzgrconsirtfit.pml(d)), error=function(e) e, warning=function(w) w)
    }

    ## Requires estfun.ltm/estfun.ltm
    if(!is(mz, "warning") & !is(mz, "try-error")){
      ord_gefp <-vector("list",length(parnum))
      for (j in 1: length(parnum)){
        if (estimator == "MML"){
        ord_gefp[[j]]  <- try(gefp(mz, fit = NULL, vcov = info_full,
                               order.by = d$age, sandwich =FALSE,
                               parm = eval(parse(text=parnum[j]))),
                          silent = TRUE)
      }
      if (estimator == "PML"){
        ord_gefp[[j]]  <- try(gefp(mz, fit = NULL, vcov = NULL,
                               order.by = d$age, sandwich =FALSE,
                               parm = eval(parse(text=parnum[j]))),
                          silent = TRUE)
      }

      if(!inherits(ord_gefp[[j]], "try-error")) {
        pval[i, (j-1)*3+1] <- sctest(ord_gefp[[j]], functional = ordfun)$p.value
        pval[i, (j-1)*3+2] <- sctest(ord_gefp[[j]], functional = ordwmax(ord_gefp[[j]]))$p.value
        pval[i, (j-1)*3+3] <- sctest(ord_gefp[[j]], functional = catL2BB(ord_gefp[[j]]))$p.value
      }
    }
   }
 }
  rval <- colMeans(pval < size, na.rm = TRUE)
  if(verbose) print(rval)

  return(rval)
}

## Loop over scenarios
simulation <- function(diff = seq(0, 4, by = 0.25),
                       parms=c("alpha", "gamma"),
                       nobs = c(120, 480, 960), nlevels = c(4, 8, 12), estimator = "PML",
                       meanch = 1, meantest = 1, verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs, parms = parms)
  nprs <- nrow(prs)
  if (meantest==1|meantest==2|meantest==0){
    parnum <- c(1,2,3,4,5,6,7,8,9,10,"1:5","6:10")
  } else{
    parnum <- c("1:5", "6:10", 11, 12)
  }

  tname <- c("ordmax","ordwmax","catdiff") # had suplm here

  test <- paste(rep(tname,length(parnum)), rep(parnum,
                                        each=length(tname)),sep="")
  ntest <- length(test)

  ## Only simulate critical values if we need it.
  if (file.exists("critvals.rda")){
    load("critvals.rda")
    ## Check to make sure it has what we need; if not,
    ## rerun critvals.
    if (!all(nlevels %in% as.numeric(names(cval)))) critvals(nlevels)
  } else {
    cval <- critvals(nlevels)
    save(cval, file="critvals.rda")
  }
  cval.conds <- as.numeric(names(cval))

  pow <- matrix(rep(NA,ntest*nprs),ncol=ntest)
  do.parallel <- require("parallel")
#  do.parallel=FALSE
  if (do.parallel){
    pow <- mclapply(1:nprs, function(i){
      testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                parms = prs$parms[i],
                parnum = parnum, nlevels = prs$nlevels[i],
                test = test, estimator = estimator, meanch = meanch,
                meantest = meantest,
                ordfun = cval[[which(cval.conds == prs$nlevels[i])]],
                verbose = verbose, ...)},
                    mc.cores = round(.75*detectCores()),
                    mc.preschedule = FALSE)
    pow <- t(matrix(unlist(pow), ntest, nprs))
  } else {
    pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    for(i in 1:nprs) {
      if(verbose) print(prs[i,])

      ## Find critical values we need
      ordfun <- cval[[which(cval.conds == prs$nlevels[i])]]

      pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                           parms = prs$parms[i],parnum = parnum,
                           nlevels = prs$nlevels[i], ordfun = ordfun,
                           test = test, estimator = estimator, meanch = meanch,
                           meantest = meantest, verbose = verbose, ...)
    }
  }

  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(rep(tname,length(parnum)), each = nprs),
                      levels=c("ordmax","ordwmax","catdiff"))
  rval$pars <- factor(rep(rep(parnum,each=length(tname)),each=nprs),
                      levels=parnum)
  rval$nobs <- factor(rval$nobs)
  rval$parms <- factor(rval$parms)
  rval$power <- as.vector(pow)
  return(rval)
}

## Simulate critical values
critvals <- function(nlevels = c(4, 8, 12), verbose = TRUE, ...)
{
  do.parallel <- require("parallel")
  if (do.parallel){
    cval <- mclapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000),
                     mc.cores = length(nlevels))
  } else {
    cval <- lapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000))
  }
  names(cval) <- nlevels

  save(cval, file="critvals.rda")
}



## Examples of running the simulation
if(FALSE){
library("strucchange")
library("mvtnorm")
source("simpara.R")

## ltm
library("ltm")
library("mirt")
source("mz.mml.R")

## pml
install.packages('lavaan_0.5-17.tar.gz', repos=NULL)
library("lavaan")
source("mz.pml.R")

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## pml
## data test
d <- dgp1(nobs = 2400, diff = 3, nlevels = 8, parms="alphaall", meanch=2)
mz <- mzgrconsirtfit.pml(d)
ord_gefp <- gefp(mz, fit = NULL, vcov = NULL,
                               order.by = d$age, sandwich = FALSE, parm =12)

pval <- sctest(ord_gefp, functional = ordwmax(ord_gefp))
sctest(mz, order.by=d$age, functional = "WDMo", plot=TRUE, parm=12)

## mml
## data test
d <- dgp1(nobs = 2400, diff = 3, nlevels = 8, parms="alphaall", meanch=2)
mz <- mzirtfit.mirt(d)
ord_gefp <- gefp(mz, fit = NULL, scores = estfun.mirt, vcov = NULL,
                               order.by = d$age, sandwich = FALSE, parm =12)
pval <- sctest(ord_gefp, functional = ordwmax(ord_gefp))
sctest(mz, order.by=d$age, functional = "WDMo", scores=estfun.mirt, vcov = NULL,
       plot=TRUE, parm=12)

## Simulation test

sim1 <- simulation(sim = "sim1", estimator = "PML", nobs = c(1200, 4800, 9600),
                   nrep = 30,
                   diff = seq(0,4,1),
                   nlevels = 8,
                   parms = c("mean"), meanch = 3,
                   meantest = 3, size=0.05)
save(sim1, file="sim4.rda")

## Short versions:
## Simulation 1:
sim1 <- simulation(sim = "sim1", estimator = "PML", nobs = c(9600), nrep = 30,
                   diff = c(4),
                   nlevels = c(8),
                   parms = c("gamma"), meanch=0, meantest=2)
save(sim1, file="sim1.rda")

sim1$test <- factor(as.character(sim1$test),
                    levels = c("ordmax","ordwmax","catdiff"),
                    labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(gamma[3]), expression(alpha[3]),
             expression(gamma[2]), expression(alpha[2]),
             expression(list(gamma[1], ldots, gamma[5])),
             expression(list(alpha[1], ldots, alpha[5])))
levels(sim1$pars) <- c("Gamma3", "Alpha3", "Gamma2", "Alpha2",
                       "All Gammas", "All Alphas")
 sim1$nlevels <- factor(sim1$nlevels)
  levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")

  sim1$nobs <- factor(sim1$nobs)
  levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")


library("lattice")
png("gamma.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "gamma" &
		 pars %in% c("Gamma3", "Alpha3",
                             "All Alphas", "All Gammas") &
                 diff %in% c(1,2,3,4) & nobs =="n=960"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", gamma[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(1, 2, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

png("mu.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | pars, group = ~ test, data = sim1,
       subset = (parms == "mean" & nobs == 4800),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", mu[4-8], ")")),
       ylab = "Power", key = mykey, as.table = TRUE)
dev.off()


library("lattice")
png("gammanobs.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
       subset = (parms == "gamma" &
		 pars %in% c("Gamma3", "Gamma2",
                             "All Gammas", "All Alphas") &
                 diff %in% c(1,2,3,4) & test =="WDM_o"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", gamma[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(1, 3, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

png("alpha.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "alpha" &
		 pars %in% c("Alpha3", "Alpha2",
                             "All Gammas", "All Alphas") &
                 diff %in% c(1,2,3,4) & nobs == "n=960"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", alpha[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(2, 4, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

library("lattice")
png("alphanobs.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
       subset = (parms == "alpha" &
		 pars %in% c("Alpha3", "Alpha2",
                             "All Gammas", "All Alphas") &
                 diff %in% c(1,2,3,4) & test =="WDM_o"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", alpha[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(2, 4, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()



## PML results organization. Parameter order is different.

sim1$test <- factor(as.character(sim1$test),
                    levels = c("ordmax","ordwmax","catdiff"),
                    labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(alpha[3]), expression(gamma[3]),
             expression(alpha[2]), expression(gamma[2]),
             expression(list(alpha[1], ldots, alpha[5])),
             expression(list(gamma[1], ldots, gamma[5])))
levels(sim1$pars) <- c("Alpha3", "Gamma3", "Alpha2", "Gamma2",
                       "All Alphas", "All Gammas")
 sim1$nlevels <- factor(sim1$nlevels)
  levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")

  sim1$nobs <- factor(sim1$nobs)
  levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")

library("lattice")
png("gamma.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "gamma" &
		 pars %in% c("Gamma3", "Gamma2",
                             "All Gammas", "All Alphas") &
                 diff %in% c(1,2,3,4) & nobs =="n=960"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", gamma[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(2, 4, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

library("lattice")
png("gammanobs.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
       subset = (parms == "gamma" &
		 pars %in% c("Gamma3", "Gamma2",
                             "All Alphas", "All Gammas") &
                 diff %in% c(1,2,3,4) & test =="WDM_o"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", gamma[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(2, 4, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

png("alpha.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "alpha" &
		 pars %in% c("Alpha3", "Alpha2",
                             "All Alphas", "All Gammas") &
                 diff %in% c(1,2,3,4) & nobs == "n=960"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", alpha[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(1, 3, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()

library("lattice")
png("alphanobs.png")
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
       subset = (parms == "alpha" &
		 pars %in% c("Alpha3", "Alpha2",
                             "All Alphas", "All Gammas") &
                 diff %in% c(1,2,3,4) & test =="WDM_o"),
       type = "b",
       xlab = expression(paste("Violation Magnitude(", alpha[3], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels =
                           parlabs[c(1, 3, 5, 6)], ...)
	   } else {
	       strip.default(which.given, factor.levels =
                             factor.levels, ...)
	   }
       })
dev.off()
}
