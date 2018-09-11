###################
## Preliminaries ##
###################

## packages
library("lavaan")
library("strucchange")
library("mvtnorm")
library("lattice")

## auxiliary code
source("mz-ordinal.R")
source("estfun-lavaan.R")
source("efpFunctional-cat.R")
source("sim-ordinal.R")

## convenience function for plotting boundaries without color
get_boundary <- function(obj, fun, trim = 0) {
  bound <- fun$boundary(0:obj$nobs/obj$nobs)
  bound <- fun$computeCritval(0.05, ncol(obj$process)) * bound
  bound <- zoo(bound, time(obj$process))
  if(trim > 0) bound <- head(tail(bound, - floor(trim * obj$nobs)), - floor(trim * obj$nobs))
  return(bound)
}

## convenience function for computing information matrix
info_full <- function(x, ...) solve(vcov(x) * nobs(x))

## for lavaan version 0.5-18 or later
info_full <- function(x, ...) solve(vcov(x, remove.duplicated = TRUE) * nobs(x))


##################
## Simulation 1 ##
##################

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1163)

## run simulation
sim1 <- simulation()

## change labels for plotting
sim1$nlevels <- factor(sim1$nlevels)
levels(sim1$nlevels) <- paste("m=",levels(sim1$nlevels),sep="")
sim1$nobs <- factor(sim1$nobs)
levels(sim1$nobs) <- paste("n=",levels(sim1$nobs),sep="")
levels(sim1$test) <- c("LM_uo", "LRT", "LM_o", "WDM_o")
sim1$test <- factor(as.character(sim1$test),
  levels = c("LM_o", "WDM_o", "LRT", "LM_uo"),
  labels = c("maxLM_o", "WDM_o", "LRT", "LM_uo"))

## visualization
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + nobs, group = ~ test, data = sim1, 
  type = "b", xlab = "Violation Magnitude", ylab = "Power", key = mykey)


##################
## Simulation 2 ##
##################

## run simulation
sim2 <- simulation(diff = seq(0, 3, by = 0.5),
  nobs = c(1200, 4800, 9600), anomaly = TRUE)

## change labels for plotting
sim2$nlevels <- factor(sim2$nlevels)
levels(sim2$nlevels) <- paste("m=",levels(sim2$nlevels),sep="")
sim2$nobs <- factor(sim2$nobs)
levels(sim2$nobs) <- paste("n=",levels(sim2$nobs),sep="")
levels(sim2$test) <- c("LM_uo", "LRT", "LM_o", "WDM_o")
sim2$test <- factor(as.character(sim2$test),
  levels = c("LM_o", "WDM_o", "LRT", "LM_uo"),
  labels = c("maxLM_o", "WDM_o", "LRT", "LM_uo"))

## visualization
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + nobs, group = ~ test, data = sim2, 
  type = "b", xlab = "Violation Magnitude", ylab = "Power", key = mykey)


##################################################
## Youth Gratitude Data from Froh et al. (2011) ##
##################################################

data("YouthGratitude", package = "psychotools")
## remove cases with 'imputed' values (not in 1, ..., 9)
yg <- YouthGratitude[apply(YouthGratitude[, 4:28], 1, function(x) all(x %in% 1:9)), ]

###################
## GQ-6 Modeling ##
###################

## tau-equivalent model
m1 <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, group = "agegroup", meanstructure=TRUE,
  group.equal = "loadings")
## parallel model
m1.f <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, group = "agegroup", meanstructure = TRUE)
## likelihood ratio test
anova(m1.f, m1)

## empirical fluctuation process for factor loadings in tau-equivalent model
gefp1 <- gefp(m1, fit = NULL, order.by = as.numeric(yg$agegroup), 
  vcov = info_full, sandwich = FALSE, parm = 1:4)

## measurement invariance tests
set.seed(1163)
ol2bb1.fun <- ordL2BB(gefp1)
sctest(gefp1, functional = ol2bb1.fun)
sctest(gefp1, functional = ordwmax(gefp1))

## visualizations
plot(gefp1, functional = ordwmax(gefp1), axes = FALSE, main = "",
     xlab = "Age group", ylab = "Weighted max statistics")
axis(1, at = 1:5, labels = c("10-11", "12-13", "14", "15", "16"))
axis(2)

plot(gefp1, functional = ol2bb1.fun, axes = FALSE, ylim = c(0, 14), main = "",
     xlab = "Age group", ylab = "LM statistics")
axis(1, at = 1:5, labels = c("10-11", "12-13", "14", "15", "16"))
axis(2)


##################
## GAC modeling ##
##################

## parallel model
m2 <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3
   f1 ~ 0*1',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = c("loadings", "residuals", "lv.variances"))
## tau-equivalent model
m2.f <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3
   f1 ~ 0*1',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")
## likelihood ratio test
anova(m2.f, m2)

## empirical fluctuation process for variances in parallel model
gefp2 <- gefp(m2, fit = NULL, order.by = as.numeric(yg$agegroup), 
  vcov = info_full, sandwich = FALSE, parm = 3:6)

## measurement invariance tests
set.seed(1163)
ol2bb2.fun <- ordL2BB(gefp2)
sctest(gefp2, functional = ol2bb2.fun)
sctest(gefp2, functional = ordwmax(gefp2))

## visualizations
plot(gefp2, functional = ordwmax(gefp2), axes = FALSE, main = "",
     xlab = "Age group", ylab = "Weighted max statistics")
axis(1, at = 1:5, labels = c("10-11", "12-13", "14", "15", "16"))
axis(2)

plot(gefp2, functional = ol2bb2.fun, axes = FALSE, ylim = c(0, 115), main = "",
     xlab = "Age group", ylab = "LM statistics")
axis(1, at = 1:5, labels = c("10-11", "12-13", "14", "15", "16"))
axis(2)

## "intermediate" model
m2.2 <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3
   f1 ~ 0*1
   gac_1 ~~ c("v1.1","v1.1","v1.2","v1.2","v1.2","v1.2")*gac_1
   gac_2 ~~ c("v2.1","v2.1","v2.2","v2.2","v2.2","v2.2")*gac_2
   gac_3 ~~ c("v3.1","v3.1","v3.2","v3.2","v3.2","v3.2")*gac_3
   f1 ~~ c("v4.1","v4.1","v4.2","v4.2","v4.2","v4.2")*f1', 
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = c("loadings", "residuals", "lv.variances"))

## likelihood ratio tests
anova(m2.2, m2)
anova(m2.f, m2.2)

## measurement invariance tests
gefp2.2 <- gefp(m2.2, fit = NULL, order.by = as.numeric(yg$agegroup), 
  vcov = info_full, sandwich = FALSE, parm = c(3:6, 13:16))

set.seed(1163)
ol2bb22.fun <- ordL2BB(gefp2.2)
sctest(gefp2.2, functional = ol2bb22.fun)
sctest(gefp2.2, functional = ordwmax(gefp2.2))


