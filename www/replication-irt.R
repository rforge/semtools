###################################################
### Preliminaries
###################################################

## load lavaan 0.5-17:
if(file.exists("~/old_package_versions/lavaan")){
  library("lavaan", lib.loc="~/old_package_versions/")
} else {
  library("lavaan")
}

## other packages
library("strucchange")
library("mvtnorm")
library("lattice")
library("mirt")

## code
source("simpara.R")
source("mz.pml.R")

## data for application
load("irtdata.rda")

###################################################
### Simulation 1
###################################################

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
if(file.exists("sim2.rda")) {
    load("sim2.rda")
} else {
    sim1 <- simulation(sim = "sim1", nobs = c(120,480,960), nrep = 5000,
                       diff = seq(0,4,1), nlevels = c(4, 8, 12),
                       parms = c("alpha", "gamma"), estimator = "PML",
                       meanch=0, meantest=0)
    save(sim1, file="sim2.rda")
}

## arrange results for display
sim1$test <- factor(as.character(sim1$test),
levels = c("ordmax","ordwmax","catdiff"),
labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(alpha[1]), expression(alpha[2]),
expression(alpha[3]), expression(alpha[4]),expression(alpha[5]),
expression(gamma[1]), expression(gamma[2]), expression(gamma[3]),
expression(gamma[4]), expression(gamma[5]),
expression(list(alpha[1], ldots, alpha[5])),
expression(list(gamma[1], ldots, gamma[5])))
levels(sim1$pars) <- c("Alpha1", "Alpha2", "Alpha3", "Alpha4", "Alpha5",
                       "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5",
"All Alphas", "All Gammas")
sim1$nlevels <- factor(sim1$nlevels)
levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")
sim1$nobs <- factor(sim1$nobs)
levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")


###################################################
### Figure 1
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & nobs == "n=960"),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 3, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 2
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha2", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & nobs =="n=960"),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(3, 7, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 3
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & test =="WDM_o"),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 3, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 4
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha3", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & test =="WDM_o"),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(3, 7, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Simulation 2.1
###################################################
## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
if(file.exists("sim8.rda")) {
load("sim8.rda")
} else {
    sim1 <- simulation(sim = "sim1", nobs = c(1200, 4800, 9600), nrep = 5000,
                       diff = seq(0,4,1), nlevels = 8,
                       parms = c("alpha", "gamma"), estimator = "PML",
                       meanch=0, meantest=2)
    save(sim1, file="sim8.rda")
}

## arrange results for display
sim1$test <- factor(as.character(sim1$test),
levels = c("ordmax","ordwmax","catdiff"),
labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(alpha[1]), expression(alpha[2]),
expression(alpha[3]), expression(alpha[4]),expression(alpha[5]),
expression(gamma[1]), expression(gamma[2]), expression(gamma[3]),
expression(gamma[4]), expression(gamma[5]),
expression(list(alpha[1], ldots, alpha[5])),
expression(list(gamma[1], ldots, gamma[5])))
levels(sim1$pars) <- c("Alpha1", "Alpha2", "Alpha3", "Alpha4", "Alpha5",
                       "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5",
"All Alphas", "All Gammas")
sim1$nlevels <- factor(sim1$nlevels)
levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")
sim1$nobs <- factor(sim1$nobs)
levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")


###################################################
### Figure 5
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nobs + pars, group = ~ test, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4)),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 3, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 6
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nobs + pars, group = ~ test, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha3", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4)),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(3, 7, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Simulation 2.2
###################################################
## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
if(file.exists("sim4.rda")) {
    load("sim4.rda")
} else {
    sim1 <- simulation(sim = "sim1", nobs = c(1200, 4800, 9600), nrep = 5000,
                       diff = seq(0,4,1), nlevels = 8,
                       parms = c("alpha", "gamma"), estimator = "PML",
                       meanch=2, meantest=2)
    save(sim1, file="sim4.rda")
}

## arrange results for display
sim1$test <- factor(as.character(sim1$test),
levels = c("ordmax","ordwmax","catdiff"),
labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(alpha[1]), expression(alpha[2]),
expression(alpha[3]), expression(alpha[4]),expression(alpha[5]),
expression(gamma[1]), expression(gamma[2]), expression(gamma[3]),
expression(gamma[4]), expression(gamma[5]),
expression(list(alpha[1], ldots, alpha[5])),
expression(list(gamma[1], ldots, gamma[5])))
levels(sim1$pars) <- c("Alpha1", "Alpha2", "Alpha3", "Alpha4", "Alpha5",
                       "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5",
"All Alphas", "All Gammas")
sim1$nlevels <- factor(sim1$nlevels)
levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")
sim1$nobs <- factor(sim1$nobs)
levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")


###################################################
### Figure 7
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nobs + pars, group = ~ test, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4)),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 3, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 8
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nobs + pars, group = ~ test, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha3", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4)),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(3, 7, 8, 11, 12)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Arrange data for application section
###################################################
ISIlevel <- unique(datstu$ISI)
SESlevel <- unique(datstu$SES)
datstu$total <- apply(datstu[,1:18], 1, mean)
datstu$SESgroup <- rep(NA, nrow = datstu)
datstu$SESgroup[datstu$SES <= -1] <- 1
datstu$SESgroup[datstu$SES > -1 & datstu$SES <= -0.5] <- 2
datstu$SESgroup[datstu$SES > -0.5 & datstu$SES <= 0] <- 3
datstu$SESgroup[datstu$SES > 0 & datstu$SES <= 0.5] <- 4
datstu$SESgroup[datstu$SES > 0.5 & datstu$SES <= 1] <- 5
datstu$SESgroup[datstu$SES > 1 & datstu$SES <= 3] <- 6
datstu$SESgroup <- as.factor(datstu$SESgroup)
datstu <- datstu[order(datstu$SESgroup),]
data <- datstu


###################################################
### Compute various test statistics
###################################################
## obtain critical values
if(file.exists("maxLMo.rda")){
  load("maxLMo.rda")
} else {
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  maxLMo <- ordL2BB(data$SESgroup, nproc = 18)
  save(maxLMo, file = "maxLMo.rda")
}
if(file.exists("maxLMo1.rda")){
  load("maxLMo1.rda")
} else {
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  maxLMo1 <- ordL2BB(data$SESgroup, nproc = 1)
  save(maxLMo1, file = "maxLMo1.rda")
}
if(file.exists("maxLMo2.rda")){
  load("maxLMo2.rda")
} else {
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  maxLMo2 <- ordL2BB(data$SESgroup, nproc = 2)
  save(maxLMo2, file = "maxLMo2.rda")
}
if(file.exists("maxLMo5.rda")){
  load("maxLMo5.rda")
} else {
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  maxLMo5 <- ordL2BB(data$SESgroup, nproc = 5)
  save(maxLMo5, file = "maxLMo5.rda")
}

## define and fit models
f05 <- paste(paste("lam", 1:18, "*", names(data)[1:18],
sep = ""), collapse = " + ")
f05 <- paste("theta =~ ", f05, sep = "")
f06 <- paste(names(data)[1:18], " | th", 1:18, "*t1",
sep = "")
f06 <- paste(f06, collapse = " \n ")
f07 <- paste("alpha", 1:18, " := (lam",1:18, ")/sqrt(1-lam",
1:18, "^2)", sep = "")
f07 <- paste(f07, collapse = " \n ")
f08 <- paste("beta", 1:18, " := (-th", 1:18, ")/sqrt(1-lam",
1:18, "^2)", sep = "")
f08 <- paste(f08, collapse = " \n ")
f09 <- paste(names(data)[1:18], " ~*~ 1*", names(data)[1:18],
               sep = "", collapse = " \n ")
f10 <- "theta ~ c(0, NA, NA, NA, NA, NA)*1"
f11 <- "theta ~~ c(1, NA, NA, NA, NA, NA)*theta"
f12 <- "theta ~ c(0, NA, NA, NA, NA, NA)*1 + c(NA, 'm', 'm', 'm', 'm', 'm')*1"
f13 <- "theta ~~ c(1, NA, NA, NA, NA, NA)*theta + c(NA, 'v', 'v', 'v', 'v', 'v')*theta"
f14 <- "theta ~~ c(1, NA, NA, NA, NA, NA)*theta + c(NA, 'v1', 'v1', 'v1', 'v2', 'v2')*theta"
mod <- paste(f05, f06, f07, f08, sep = " \n ")
writeLines(mod, "mod.txt")

modgroup <- paste(f05, f06, f07, f08, f09, f10, f11, sep = " \n ")
writeLines(modgroup, "modgroup.txt")

modmean <-  paste(f05, f06, f07, f08, f09, f10, sep = " \n ")
writeLines(modmean, "modmean.txt")

modgrcons <- paste(f05, f06, f07, f08, f09, f12, f13, sep = " \n ")
writeLines(modgrcons, "modgrcons.txt")

modgrconsfinal <- paste(f05, f06, f07, f08, f09, f10, f14, sep = " \n ")
writeLines(modgrconsfinal, "modgrconsfinal.txt")

if(file.exists("rval.rda")){
  load("rval.rda")
} else {
  rval <- cfa(mod, data=data, std.lv=TRUE,
            ordered=names(data)[1:18], estimator="PML")
  save(rval, file="rval.rda")
}
if(file.exists("rvalgroup.rda")){
  load("rvalgroup.rda")
} else {
  rvalgroup <- cfa(modgroup, data=data, std.lv=TRUE,
                 ordered=names(data)[1:18], group = "SESgroup", group.equal=
                 c("loadings", "thresholds"), estimator="PML")
  save(rvalgroup, file="rvalgroup.rda")
}
if(file.exists("rvalmean.rda")){
  load("rvalmean.rda")
} else {
  rvalmean <- cfa(modmean, data=data, std.lv=TRUE,
                 ordered=names(data)[1:18], group = "SESgroup", group.equal=
                 c("loadings", "thresholds"), estimator="PML")
  save(rvalmean, file="rvalmean.rda")
}
if(file.exists("rvalgrcons.rda")){
  load("rvalgrcons.rda")
} else {
  rvalgrcons <- cfa(modgrcons, data=data, std.lv=TRUE,
                 ordered=names(data)[1:18], group = "SESgroup", group.equal=
                 c("loadings", "thresholds"), estimator="PML")
  save(rvalgrcons, file="rvalgrcons.rda")
}
if(file.exists("rvalgrconsfinal.rda")){
  load("rvalgrconsfinal.rda")
} else {
  rvalgrconsfinal <- cfa(modgrconsfinal, data=data, std.lv=TRUE,
                 ordered=names(data)[1:18], group = "SESgroup", group.equal=
                 c("loadings", "thresholds"), estimator="PML")
  save(rvalgrconsfinal, file="rvalgrconsfinal.rda")
}
WDMoSESdisc <- rep(NA, 18)
WDMoSESdiff <- rep(NA, 18)
WDMoSESgroupdisc <- rep(NA, 18)
WDMoSESgroupdiff <- rep(NA, 18)
maxLMoSESdisc <- rep(NA, 18)
maxLMoSESdiff <- rep(NA, 18)
maxLMoSESgroupdisc <- rep(NA, 18)
maxLMoSESgroupdiff <- rep(NA, 18)

## obtain test statistics/p-values
if(file.exists("parares.rda")){
  load("parares.rda")
} else {
  for (i in 1:18){
    WDMoSESdisc[i] <- sctest(rval, order.by = data$SESgroup, parm = i,
                             vcov = NULL,
                             functional = "WDMo", plot = FALSE)$p.value
    WDMoSESdiff[i] <- sctest(rval, order.by = as.factor(data$SESgroup),
                             parm = i+18, vcov = NULL,
                             functional = "WDMo", plot = FALSE)$p.value
    WDMoSESgroupdisc[i] <- sctest(rvalmean, order.by = data$SESgroup,
                                  parm = i, vcov = NULL,
                                  functional = "WDMo", plot = FALSE)$p.value
    WDMoSESgroupdiff[i] <- sctest(rvalmean, order.by = data$SESgroup,
                                  parm = i+18, vcov = NULL,
                                  functional = "WDMo", plot = FALSE)$p.value
    maxLMoSESdisc[i] <- sctest(rval, order.by = data$SESgroup, parm = i,
                             vcov = NULL,
                             functional = maxLMo1, plot = FALSE)$p.value
    maxLMoSESdiff[i] <- sctest(rval, order.by = as.factor(data$SESgroup),
                             parm = i+18, vcov = NULL,
                             functional = maxLMo1, plot = FALSE)$p.value
    maxLMoSESgroupdisc[i] <- sctest(rvalmean, order.by = data$SESgroup,
                                  parm = i, vcov = NULL,
                                  functional = maxLMo1, plot = FALSE)$p.value
    maxLMoSESgroupdiff[i] <- sctest(rvalmean, order.by = data$SESgroup,
                                  parm = i+18, vcov = NULL,
                                  functional = maxLMo1, plot = FALSE)$p.value
  }

  parares <- cbind(WDMoSESdisc, WDMoSESdiff, WDMoSESgroupdisc,
                   WDMoSESgroupdiff, maxLMoSESdisc, maxLMoSESdiff,
                   maxLMoSESgroupdisc, maxLMoSESgroupdiff)
  save(parares, file ="parares.rda")
}

if(file.exists("personparares.rda")){
  load("personparares.rda")
} else {
    personparares <- rep(NA, 3)
    personparares[1] <- sctest(rvalgrcons, order.by = data$SESgroup,
                                  parm = 37, vcov = NULL,
                            functional = maxLMo1, plot = TRUE)$p.value
    personparares[2] <- sctest(rvalgrcons, order.by = data$SESgroup,
                                  parm = 38, vcov = NULL,
                               functional = maxLMo1, plot = FALSE)$p.value
    personparares[3] <- sctest(rvalgrconsfinal, order.by = data$SESgroup,
                                  parm = 38, vcov = NULL,
                               functional = maxLMo1, plot = FALSE)$p.value
    save(personparares, file ="personparares.rda")
}


###################################################
### Appendix B: Results for MML
###################################################
## ltm
library("ltm")
library("strucchange")
library("mvtnorm")
library("lattice")
source("mz.mml.R")

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1090)

## run simulation
if(file.exists("sim1.rda")) {
  load("sim1.rda")
} else {
sim1 <- simulation(sim = "sim1", nobs = c(120,480,960), nrep = 5000,
                       diff = seq(0,4,1), nlevels = c(4, 8, 12),
                       parms = c("alpha", "gamma"), estimator = "MML",
                       meanch=0, meantest=0)
save(sim1, file="sim1.rda")
}
sim1$test <- factor(as.character(sim1$test),
levels = c("ordmax","ordwmax","catdiff"),
labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(alpha[1]), expression(alpha[2]),
expression(alpha[3]), expression(alpha[4]),expression(alpha[5]),
expression(gamma[1]), expression(gamma[2]), expression(gamma[3]),
expression(gamma[4]), expression(gamma[5]),
expression(list(alpha[1], ldots, alpha[5])),
expression(list(gamma[1], ldots, gamma[5])))
levels(sim1$pars) <- c("Alpha1", "Alpha2", "Alpha3", "Alpha4", "Alpha5",
                       "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5",
"All Alphas", "All Gammas")
sim1$nlevels <- factor(sim1$nlevels)
levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")
sim1$nobs <- factor(sim1$nobs)
levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")


###################################################
### Figure 11
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & nobs == "n=960"),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 8, 3, 12, 11)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 12
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha3", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & nobs =="n=960"),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(8, 3, 7, 12, 11)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 13
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
subset = (parms == "alpha" &
pars %in% c("Alpha2", "Alpha3", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & test =="WDM_o"),
type = "b",
xlab = expression(paste("Violation Magnitude (", alpha[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(2, 8, 3, 12, 11)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})


###################################################
### Figure 14
###################################################
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$nobs), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ nobs, data = sim1,
subset = (parms == "gamma" &
pars %in% c("Alpha3", "Gamma2", "Gamma3",
"All Alphas", "All Gammas") &
diff %in% c(0,1,2,3,4) & test =="WDM_o"),
type = "b",
xlab = expression(paste("Violation Magnitude (", gamma[3], ")")),
ylab = "Power", key = mykey, as.table = TRUE,
strip = function(..., which.given, factor.levels){
if(which.given == 2){
strip.default(which.given, factor.levels =
parlabs[c(8, 3, 7, 12, 11)], ...)
} else {
strip.default(which.given, factor.levels =
factor.levels, ...)
}
})
