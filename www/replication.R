###################
## Preliminaries ##
###################

## packages and code
library("OpenMx")
library("lavaan")
library("strucchange")
library("mvtnorm")
library("lattice")
source("mz.R")
source("estfun-lavaan.R")
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
  vcov = info.mzfit, sandwich = FALSE)
gefp_19_info <- gefp(mz0, fit = NULL, order.by = d$age,
  vcov = info.mzfit, sandwich = FALSE)

## Figure 2
par(mfcol = c(3, 2))
plot(gefp_3_info,  functional = maxBB, main = "DM, k* = 3", xlab = "Age")
plot(gefp_3_info,  functional = meanL2BB, main = "CvM, k* = 3", xlab = "Age")
plot(gefp_3_info,  functional = supLM(0.1), main = "max LM, k* = 3", xlab = "Age")
plot(gefp_19_info, functional = maxBB, main = "DM, k* = 19", xlab = "Age")
plot(gefp_19_info, functional = meanL2BB, main = "CvM, k* = 19", xlab = "Age")
plot(gefp_19_info, functional = supLM(0.1), main = "max LM, k* = 19", xlab = "Age")


## Figure 3
## compute maxLR "by hand"
lrstat <- function(age) 2 * as.vector(logLik(mzfit(d[d$age <= age,])) + logLik(mzfit(d[d$age > age,])) - logLik(mz0))
mz.age <- sort(d$age[d$age >= quantile(d$age, 0.1) & d$age <= quantile(d$age, 0.9)])
mz.lrstat <- sapply(mz.age, lrstat)

par(mfrow = c(1, 1))
plot(gefp_19_info, functional = supLM(0.1),
  xlab = "Age", ylab = "LR and LM statistics (k* = 19)", main = "", ylim = c(0, 50))
legend("topleft", c("LR", "LM"), lty = c(2, 1), bty = "n")
lines(mz.age, mz.lrstat, lty = 2)
abline(v = mz.age[which.max(mz.lrstat)], lty = 3)


## Figure 4
plot(gefp_3_info,  functional = maxBB, aggregate = FALSE, ylim = c(-2, 2),
  main = "DM, k* = 3", xlab = "Age", ylab = expression(lambda[11], lambda[21], lambda[31]))


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

## display subset of results in table
## (with rounding, in percent)
mz_sim_sub <- subset(mz_sim, pars %in% c(3, 19))
mz_sim_sub$pars <- factor(mz_sim_sub$pars)
levels(mz_sim_sub$pars) <- paste("k* =", levels(mz_sim_sub$pars))
levels(mz_sim_sub$nobs) <- paste("n =", levels(mz_sim_sub$nobs))
round(ftable(100 * xtabs(power ~ nobs + pars + test + diff,
  data = mz_sim_sub, subset = diff %in% c(seq(0, 3.5, by = 0.5))),
  col.vars = "diff"), digits = 1)

## display result in graphic
library("lattice")
trellis.par.set(theme = canonical.theme(color = FALSE))
xyplot(power ~ diff | pars + nobs, group = ~ test, data = mz_sim_sub, type = "b",
  xlab = "Violation Magnitude", ylab = "Power", ylim = c(0, 1))


####################################
## Application: Stereotype Threat ##
####################################

## data and transformations
data("StereotypeThreat", package = "psychotools")
## include group variable
StereotypeThreat <- transform(StereotypeThreat, group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]
## omit NAs in GPA and break ties randomly
StereotypeThreat <- subset(StereotypeThreat, !is.na(gpa))
ties <- duplicated(StereotypeThreat$gpa)
set.seed(1090)
StereotypeThreat$gpa[ties] <- runif(sum(ties),
  min = StereotypeThreat$gpa[ties] - 0.001, max = StereotypeThreat$gpa[ties] + 0.01)
## keep only data without NAs in GPA
wdh <- subset(StereotypeThreat, !is.na(gpa))


## fit model 5b
wdh5b <- lavaan(
   'ability   =~ label(c(rep("load_n", 3), "load_n:min_t")) * numerical + label(rep("load_v", 4)) * verbal + 1 * abstract
    ability   ~  label(c(NA, "lmean:min_c", "lmean:maj_t", "lmean:min_t")) * 1 + c(0, NA, NA, NA) * 1
    ability   ~~ label(c("lvar:maj", "lvar:min", "lvar:maj", "lvar:min")) * ability
    numerical ~  label(c(rep("mean_n", 3), "mean_n:min_t")) * 1
    abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", rep("mean_a", 2))) * 1
    verbal    ~  label(rep("mean_v", 4)) * 1
    numerical ~~ label(c(rep("var_n", 3), "var_n:min_t")) * numerical
    abstract  ~~ label(rep("var_a", 4)) * abstract
    verbal    ~~ label(rep("var_v", 4)) * verbal',
  data = wdh, meanstructure = TRUE, group = "group")

## refit model with new identification constraints
wdh5b_1 <- lavaan(
   'ability   =~ 1 * verbal + label(c(rep("load_n", 3), "load_n:min_t")) * numerical  + label(rep("load_a", 4)) * abstract
    ability   ~  label(c(NA, "lmean:min_c", "lmean:maj_t", "lmean:min_t")) * 1 + c(0, NA, NA, NA) * 1
    ability   ~~ label(c("lvar:maj", "lvar:min", "lvar:maj", "lvar:min")) * ability
    numerical ~  label(c(rep("mean_n", 3), "mean_n:min_t")) * 1
    abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", rep("mean_a", 2))) * 1
    verbal    ~  label(rep("mean_v", 4)) * 1
    numerical ~~ label(c(rep("var_n", 3), "var_n:min_t")) * numerical
    abstract  ~~ label(rep("var_a", 4)) * abstract
    verbal    ~~ label(rep("var_v", 4)) * verbal',
  data = wdh, meanstructure = TRUE, group = "group")


## extract information matrix and scores in the minority/threat group
info_full <- function(x, ...) solve(vcov(x) * nobs(x))
scores_min_t <- function(x, ...) {
  ef <- estfun(x, ...)
  ef[wdh$group == "minority.threat", ]
}
gpa_min_t <- subset(wdh, group == "minority.threat")$gpa

## fluctuation processes for the four group-specific parameters
## in both models (i.e., with different identification constraints)
gefp_4_info <- gefp(wdh5b, fit = NULL, scores = scores_min_t,
  order.by = gpa_min_t, vcov = info_full, sandwich = FALSE, parm = 15:18)
gefp_4_info_1 <- gefp(wdh5b_1, fit = NULL, scores = scores_min_t,
  order.by = gpa_min_t, vcov = info_full, sandwich = FALSE, parm = 15:18)


## Figure 7
plot(gefp_4_info,  functional = maxBB, aggregate = FALSE,
  ylim = c(-2.15, 2.15), main = "DM, k* = 4", xlab = "GPA",
  ylab = expression(lambda[num], eta, mu[num], psi[num]))

## p-values from all three tests under both sets of constraints
tests <- list(maxBB, meanL2BB, supLM(0.1))
pvals <- matrix(NA, nrow = length(tests), ncol = 2)
for (i in 1:length(tests)){
  pvals[i, 1] <- sctest(gefp_4_info,   functional = tests[[i]])$p.value
  pvals[i, 2] <- sctest(gefp_4_info_1, functional = tests[[i]])$p.value
}
