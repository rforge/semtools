###################
## Preliminaries ##
###################

## packages and code
library("OpenMx")
library("strucchange")
library("mvtnorm")
source("mz.R")
source("sim.R")

##################################
## Example with Artificial Data ##
##################################

## generate single artificial data set
set.seed(1090)
d <- dgp(200, 2)

## full sample and "true" subsample models 
mz0 <- mzfit(d)
mz1 <- mzfit(d[1:100,])
mz2 <- mzfit(d[101:200,])

## check: Which model parameters differ pre- and post-break?
round((coef(mz1) - coef(mz2)) / sqrt(diag(vcov(mz1)) + diag(vcov(mz2))), digits = 2)

## empirical fluctuation processes (3 vs all 19 pars)
gefp_3_info <- gefp(mz0, fit = NULL, order.by = d$age, parm = 1:3,
  vcov = info.MxModel, sandwich = FALSE)
gefp_19_info <- gefp(mz0, fit = NULL, order.by = d$age,
  vcov = info.MxModel, sandwich = FALSE)

## Figure 2
 # pdf("gefp_3_19.pdf",height = 9, width = 6)
par(mfcol = c(3, 2))
plot(gefp_3_info,  functional = maxBB, main = "Three Parameters\nDM", xlab = "Age")
plot(gefp_3_info,  functional = meanL2BB, main = "CvM", xlab = "Age")
plot(gefp_3_info,  functional = supLM(0.1), main = "maxLM", xlab = "Age")
plot(gefp_19_info, functional = maxBB, main = "Nineteen Parameters\nDM", xlab = "Age")
plot(gefp_19_info, functional = meanL2BB, main = "CvM", xlab = "Age")
plot(gefp_19_info, functional = supLM(0.1), main = "maxLM", xlab = "Age")
 # dev.off()


## Figure 3
## compute maxLR "by hand"
lrstat <- function(age) 2 * as.vector(logLik(mzfit(d[d$age <= age,])) + logLik(mzfit(d[d$age > age,])) - logLik(mz0))
mz.age <- sort(d$age[d$age >= quantile(d$age, 0.1) & d$age <= quantile(d$age, 0.9)])
mz.lrstat <- sapply(mz.age, lrstat)

 # pdf("lmlr_19.pdf", height = 5, width = 6)
par(mfrow = c(1, 1))
plot(gefp_19_info, functional = supLM(0.1),
  xlab = "Age", ylab = "LR and LM statistics", main = "", ylim = c(0, 50))
legend("topleft", c("LR", "LM"), lty = c(1, 2), bty = "n")
lines(mz.age, mz.lrstat, lty = 2)
abline(v = mz.age[which.max(mz.lrstat)], lty = 3)
 # dev.off()


## Figure 4
 # pdf("gefp_lambdas.pdf", height = 7, width = 6)
plot(gefp_3_info,  functional = maxBB, aggregate = FALSE, ylim = c(-2, 2),
  main = "", xlab = "Age", ylab = expression(lambda[11], lambda[21], lambda[31]))
 # dev.off()


## Find points at which test statistics achieve their maxima:
# L2 statistics:
L2_process <- zoo(apply(as.matrix(gefp_3_info$process)^2, 1, sum),
                    time(gefp_3_info))
L2_max <- L2_process[which.max(L2_process)]
# Double max statistic:
dmax_process <- zoo(apply(abs(as.matrix(gefp_3_info$process)), 1, max),
                    time(gefp_3_info))
dmax <- dmax_process[which.max(dmax_process)]


################
## Simulation ##
################

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
mz_sim <- simulation()

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
levels(mz_sim$pars) <- paste("k =", levels(mz_sim$pars))
levels(mz_sim$nobs) <- paste("n =", levels(mz_sim$nobs))
xyplot(power ~ diff | pars + nobs, group = ~ test, data = mz_sim, type = "b",
       xlab = "Violation Magnitude", ylab = "Power", ylim = c(0, 1))
