## packages and code
library("OpenMx")
library("strucchange")
library("mvtnorm")
source("mz.R")
source("sim.R")

## generate artificial data
set.seed(1090)
d <- dgp(200, 2)

## full sample and "true" subsample models 
mz0 <- mzfit(d)
mz1 <- mzfit(d[1:100,])
mz2 <- mzfit(d[101:200,])

## check: Which model parameters differ pre- and post-break?
(coef(mz1) - coef(mz2)) / sqrt(diag(vcov(mz1)) + diag(vcov(mz2)))

## empirical fluctuation processes (3 vs 19 pars, info vs OPG vcov)
gefp_3_opg <- gefp(mz0, fit = NULL, order.by = d$age, parm = 1:3)
gefp_19_opg <- gefp(mz0, fit = NULL, order.by = d$age)
gefp_3_info <- gefp(mz0, fit = NULL, vcov = info.MxModel,
  order.by = d$age, sandwich = FALSE, parm = 1:3)
gefp_19_info <- gefp(mz0, fit = NULL, vcov = info.MxModel,
  order.by = d$age, sandwich = FALSE)

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
lrstat <- function(age) 2 * as.vector(logLik(mzfit(d[d$age <= age,])) + logLik(mzfit(d[d$age > age,])) - logLik(mz0))
mz.age <- sort(d$age[d$age >= quantile(d$age, 0.1) & d$age <= quantile(d$age, 0.9)])
mz.lrstat <- sapply(mz.age, lrstat)

## compare supLM and supLR test
## (the critical value is the same)
par(mfrow = c(1, 1))
plot(gefp_19_info, functional = supLM(0.1),
  xlab = "Age", ylab = "LM and LR statistics", main = "",
  ylim = c(0, 65))
lines(mz.age, mz.lrstat, col = 4)
abline(v = mz.age[which.max(mz.lrstat)], lty = 2)

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


library("strucchange")
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)
sc_sim <- simulation()

## Results in tabular form (replication of Table 6)
tab <- xtabs(power ~ intensity + test + angle + timing, data = sc_sim)
ftable(tab, row.vars = c("angle", "timing", "test"), col.vars = "intensity")
