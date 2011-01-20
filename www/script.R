# Main R script file to carry out tests of measurement invariance
# described in Merkle & Zeileis.

# Load requisite packages
library("OpenMx")
library("strucchange")
library("numDeriv")

# Read data
obs.data <- dget("data.dat")
dat <- obs.data$datmat
V <- obs.data$age

# Obtain functions for score and information matrix extraction
source("fa_score_extraction.R")

# Obtain OpenMx model specification
# (Two-group model with all parameters constrained to be equal,
#  though this could also be estimated as a single-group model.)
source("def.model.R")
fa.mod <- def.model(dat, V)

# Fit model via OpenMx
mx.res <- mxRun(fa.mod)

## the same with a single model (not split/constrained by age)
fa.mod1 <- def.model1(dat)
mx.res1 <- mxRun(fa.mod1)

# Organize results for strucchange:
fa.res <- list(res = mx.res1, dat = dat, V = V)

## store scores in list as well to avoid slow re-computation
fa.res$estfun <- score.fa(fa.res)
colnames(fa.res$estfun) <- summary(mx.res1)$parameters$name

# Calculate cumulative scores for first three factor loadings:
fa.cums1 <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  vcov = inf.fa, order.by = fa.res$V, sandwich = FALSE, parm = 1:3)

## the same for all parameters
fa.cums2 <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  vcov = inf.fa, order.by = fa.res$V, sandwich = FALSE)

## the same with OPG (outer product of gradients) vcov estimator
fa.cums1opg <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  order.by = fa.res$V, parm = 1:3)
fa.cums2opg <- gefp(fa.res, fit=NULL, scores = function(x) x$estfun,
  order.by = fa.res$V)


# Plot individual cumulative score processes:
plot(fa.cums1, aggregate = FALSE, ylim = c(-1.8, 1.8), xlab = "Age")

## Cramer-von Mises test
## -> all lead to significant results but with all 19 parameters
##    the peak is even clearer at closer to age = 16
par(mfrow = c(2, 2))
plot(fa.cums1,    functional = meanL2BB, main = "CvM: 3 pars, Info")
plot(fa.cums2,    functional = meanL2BB, main = "CvM: 19 pars, Info")
plot(fa.cums1opg, functional = meanL2BB, main = "CvM: 3 pars, OPG")
plot(fa.cums2opg, functional = meanL2BB, main = "CvM: 19 pars, OPG")

## analogous results for supLM test (with 10% trimming)
par(mfrow = c(2, 2))
plot(fa.cums1,    functional = supLM(0.1), main = "supLM: 3 pars, Info")
plot(fa.cums2,    functional = supLM(0.1), main = "supLM: 19 pars, Info")
plot(fa.cums1opg, functional = supLM(0.1), main = "supLM: 3 pars, OPG")
plot(fa.cums2opg, functional = supLM(0.1), main = "supLM: 19 pars, OPG")

## compute supLR "by hand"
lrstat <- function(age) {
  summary(mx.res1)$Minus2LogLikelihood -
  summary(mxRun(def.model1(dat[V <= age,])))$Minus2LogLikelihood -
  summary(mxRun(def.model1(dat[V >  age,])))$Minus2LogLikelihood
}
fa.age <- sort(V[V >= quantile(V, 0.1) & V <= quantile(V, 0.9)])
fa.lrstat <- sapply(fa.age, lrstat)

## compare supLM and supLR test
## (the critical value is the same)
par(mfrow = c(1, 1))
plot(fa.cums2, functional = supLM(0.1),
  xlab = "Age", ylab = "LM and LR statistics", main = "",
  ylim = c(0, 65))
lines(fa.age, fa.lrstat, col = 4)
abline(v = fa.age[which.max(fa.lrstat)], lty = 2)
## -> in this case, the conclusions are exactly the same, the
##    maximum is even assumed for the same value (15.96474),
##    the remaining differences may possibly be reduced when
##    the numerical computations are improved
