####################
### Preliminaries ##
####################

## packages
library("lavaan")
library("strucchange") ## at least 1.5-0
library("mvtnorm")
library("lattice")

## auxiliary code
source("../www/sim-frontiers.R")
source("../www/mz-frontiers.R")

## some options
options(prompt = "R> ", continue = "+  ", digits = 3,
  useFancyQuotes = FALSE, show.signif.stars = FALSE)


#######################
### Tutorial Section ##
#######################

data("YouthGratitude", package = "psychotools")
compcases <- apply(YouthGratitude[, 4:28], 1, function(x) all(x %in% 1:9))
yg <- YouthGratitude[compcases, ]

restr <- cfa("f1 = ~gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5",
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")
full <- cfa("f1= ~gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5",
  data = yg, group = "agegroup", meanstructure = TRUE)

anova(full, restr)

sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = "LMuo")

dm <- sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = "DM")
cvm <- sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = "CvM")
maxlm <- sctest(restr, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = "maxLM")
c(dm$p.value, cvm$p.value, maxlm$p.value)

set.seed(1090)
wdmo <- sctest(restr, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = "WDMo")

## Use of saved critical values:
if(file.exists("maxLMo.rda")){
  load("maxLMo.rda")
} else {
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  maxLMo <- ordL2BB(yg$agegroup, nproc = 4)
  save(maxLMo, file = "maxLMo.rda")
}
maxlmo <- sctest(restr, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = maxLMo)

wdmo <- sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = "WDMo")
## Use without saving critical values:
maxlmo <- sctest(restr, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = "maxLMo")

c(wdmo$p.value, maxlmo$p.value)


#########################
### Instability graphs ##
#########################

par(mfrow = c(1, 2))
set.seed(1090)
sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = "WDMo", plot = TRUE, ylim = c(0, 3.2),
  xlab = "Age group", ylab = "Weighted test statistics",
  main = expression(WDM[o]), boundary = list(col = 1, lty = 2))
sctest(restr, order.by = yg$agegroup, parm = 1:4, vcov = "info",
  functional = maxLMo, plot = TRUE, ylim = c(0, 14),
  xlab = "Age group", ylab = "LM statistics",
  main = expression(maxLM[o]), boundary = list(col = 1, lty = 2))


############################
### Saved critical values ##
############################

mLMo <- ordL2BB(yg$agegroup)
maxlmo <- sctest(restr, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = mLMo)


###################
### Simulation 1 ##
###################

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1163)

if(file.exists("sim1.rda")) {
  load("sim1.rda")
} else {
  sim1 <- simulation(sim = "sim1", nobs = 480, parms = c("loading", "error", "var"))
  sim1$nlevels <- factor(sim1$nlevels)
  levels(sim1$nlevels) <- paste("m=", levels(sim1$nlevels), sep="")
 
  sim1$nobs <- factor(sim1$nobs)
  levels(sim1$nobs) <- paste("n=", levels(sim1$nobs), sep="")
  save(sim1, file="sim1.rda")
}
sim1$test <- factor(as.character(sim1$test), levels = c("ordmax","ordwmax","catdiff"),
                    labels = c("maxLM_o", "WDM_o", "LM_uo"))
parlabs <- c(expression(lambda[11]), expression(phi[12]), 
             expression(psi[11]), expression(list(lambda[11], ldots, lambda[62])), 
             expression(list(psi[11], ldots, psi[66])))
levels(sim1$pars) <- c("Loading1", "Covariance", "Error1", "All Loadings", "All Errors")  

## visualization
## Figure 3
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "loading" & 
		 pars %in% c("Loading1", "Error1", "All Loadings") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", lambda[11], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
           if(which.given == 2){
             strip.default(which.given, factor.levels = parlabs[c(1, 3, 4)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 4
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "var" & 
		 pars %in% c("Covariance") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", phi[12], ")")),
       ylab="Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(2)],...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 5
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim1$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim1,
       subset = (parms == "error" & 
		 pars %in% c("Error1", "All Errors") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", psi[11], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(3,5)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })


###################
### Simulation 2 ##
###################

## seed for replication
RNGkind(kind = "default", normal.kind = "default")
set.seed(1163)

if(file.exists("sim2.rda")) {
  load("sim2.rda")
} else {
  sim2 <- simulation(sim = "sim2", nobs = 480, 
  parms = c("extra", "extra+loading", "extra+var", "extra+error"))

  save(sim2, file="sim2.rda")
}
sim2$test <- factor(as.character(sim2$test),
  levels = c("ordmax", "ordwmax", "catdiff"),
  labels = c("maxLM_o", "WDM_o", "LM_uo"))
  sim2$nlevels <- factor(sim2$nlevels)
  levels(sim2$nlevels) <- paste("m=", levels(sim2$nlevels), sep="")
  sim2$nobs <- factor(sim2$nobs)
  levels(sim2$nobs) <- paste("n=", levels(sim2$nobs), sep="")
  parlabs <- c(expression(lambda[11]), expression(phi[12]), 
               expression(psi[11]),
	       expression(list(lambda[11], ldots, lambda[62])), 
               expression(list(psi[11], ldots, psi[66])))
  levels(sim2$pars) <- c("Loading1", "Covariance", "Error1", "All Loadings", "All Errors") 

## visualization
## Figure 6
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim2,
       subset = (parms == "extra" & 
		 pars %in% c("Loading1", "Covariance", "Error1") & diff %% 0.5 == 0),
       type = "b", xlab = "Violation Magnitude (Unmodeled Loading)",
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(1, 2, 3)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 7
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim2,
       subset = (parms == "extra" & 
		 pars %in% c("All Loadings", "All Errors") & diff %% 0.5 == 0),
       type = "b", xlab = "Violation Magnitude (Unmodeled Loading)", ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(4, 5)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 8
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim2,
       subset = (parms=="extra+loading" & 
		 pars %in% c("Loading1", "All Loadings") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", lambda[11], ")")), ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(1, 4)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 9
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim2,
       subset = (parms=="extra+var" & 
		 pars %in% c("Covariance") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", phi[12], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(2, 4:5)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })

## Figure 10
trellis.par.set(theme = canonical.theme(color = FALSE))
mykey <- simpleKey(levels(sim2$test), points = TRUE, lines = TRUE)
xyplot(power ~ diff | nlevels + pars, group = ~ test, data = sim2,
       subset = (parms == "extra+error" & 
		 pars %in% c("Error1", "All Errors") & diff %% 0.5 == 0),
       type = "b", xlab = expression(paste("Violation Magnitude (", psi[11], ")")),
       ylab = "Power", key = mykey, as.table = TRUE,
       strip = function(..., which.given, factor.levels){
	   if(which.given == 2){
	       strip.default(which.given, factor.levels = parlabs[c(3,5)], ...)
	   } else {
	       strip.default(which.given, factor.levels = factor.levels, ...)
	   }
       })


