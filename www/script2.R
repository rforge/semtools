## packages
library("OpenMx")
library("strucchange")
library("numDeriv")
library("mvtnorm")

## additional code
source("gen.data.R")
source("fa_score_extraction.R")
source("def.model.R")

## generate artificial data
set.seed(1090)
obs.data <- gen.data(200, 2)
dat <- obs.data$datmat
V <- obs.data$age

## full sample and "true" subsample models 
mx.res0 <- mxRun(def.model1(dat))
mx.res1 <- mxRun(def.model1(dat[1:100,]))
mx.res2 <- mxRun(def.model1(dat[101:200,]))

## check: Which model parameters differ pre- and post-break?
round(structure( 
  (summary(mx.res1)$parameters[,5] - summary(mx.res2)$parameters[,5]) /
    sqrt(summary(mx.res1)$parameters[,6]^2 + summary(mx.res2)$parameters[,6]^2),
  names = summary(mx.res1)$parameters[,1]
), digits = 2)

## organize results for strucchange
fa.res <- list(res = mx.res0, dat = dat, V = V)
fa.res$estfun <- score.fa(fa.res)
colnames(fa.res$estfun) <- summary(mx.res0)$parameters$name

## empirical fluctuation processes (3 vs 19 pars, info vs OPG vcov)
gefp_3_info <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  vcov = inf.fa, order.by = fa.res$V, sandwich = FALSE, parm = 1:3)
gefp_19_info <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  vcov = inf.fa, order.by = fa.res$V, sandwich = FALSE)
gefp_3_opg <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  order.by = fa.res$V, parm = 1:3)
gefp_19_opg <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  order.by = fa.res$V)

## Double-maximum test
par(mfrow = c(2, 2))
plot(gefp_3_info,  functional = maxBB, main = "dmax: 3 pars, Info")
plot(gefp_19_info, functional = maxBB, main = "dmax: 19 pars, Info")
plot(gefp_3_opg,   functional = maxBB, main = "dmax: 3 pars, OPG")
plot(gefp_19_opg,  functional = maxBB, main = "dmax: 19 pars, OPG")

## Cramer-von Mises test
par(mfrow = c(2, 2))
plot(gefp_3_info,  functional = meanL2BB, main = "CvM: 3 pars, Info")
plot(gefp_19_info, functional = meanL2BB, main = "CvM: 19 pars, Info")
plot(gefp_3_opg,   functional = meanL2BB, main = "CvM: 3 pars, OPG")
plot(gefp_19_opg,  functional = meanL2BB, main = "CvM: 19 pars, OPG")

## supLM test (with 10% trimming)
par(mfrow = c(2, 2))
plot(gefp_3_info,  functional = supLM(0.1), main = "supLM: 3 pars, Info")
plot(gefp_19_info, functional = supLM(0.1), main = "supLM: 19 pars, Info")
plot(gefp_3_opg,   functional = supLM(0.1), main = "supLM: 3 pars, OPG")
plot(gefp_19_opg,  functional = supLM(0.1), main = "supLM: 19 pars, OPG")

## compute supLR "by hand"
lrstat <- function(age) {
  summary(mx.res0)$Minus2LogLikelihood -
  summary(mxRun(def.model1(dat[V <= age,])))$Minus2LogLikelihood -
  summary(mxRun(def.model1(dat[V >  age,])))$Minus2LogLikelihood
}
fa.age <- sort(V[V >= quantile(V, 0.1) & V <= quantile(V, 0.9)])
fa.lrstat <- sapply(fa.age, lrstat)

## compare supLM and supLR test
## (the critical value is the same)
par(mfrow = c(1, 1))
plot(gefp_19_info, functional = supLM(0.1),
  xlab = "Age", ylab = "LM and LR statistics", main = "",
  ylim = c(0, 65))
lines(fa.age, fa.lrstat, col = 4)
abline(v = fa.age[which.max(fa.lrstat)], lty = 2)

## EMPIRICAL RESULTS (all conforming with theory!):
##   - Only first 3 of all 19 parameters change.
##   - When only the first 3 parameters are assessed, all tests
##     (dmax, CvM, supLM) have sufficient power to detect a change.
##   - When all 19 parameters are assessed, all tests have less power.
##     However, dmax is still significant, the others are not.
##   - supLR (which is by construction for all 19 parameters)
##     is close to supLM
##   - Tests with information matrix seem to be more powerful
##     than OPG.
