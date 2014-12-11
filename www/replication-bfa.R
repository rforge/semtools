
###############################################
## Preliminary packages, data, files
library("knitr")
## global knitr options
opts_chunk$set(echo=FALSE, fig.path='figure/vuong-',fig.align='center',fig.show='hold',size='footnotesize', cache.path="cache/", warning=FALSE, message = FALSE)

## packages
library("rjags", quietly=TRUE)
library("mvtnorm")
library("MCMCpack")

## auxiliary code
source("bfcalc.R")
source("jagsmodels.R")

## real-world data
load("peters08.rda")


###############################################
## Figure 1
par(mfrow=c(1,2))
shades <- paste("grey", c(100, 95, 85, 70, 60, 30))
ptshading <- shades[peters08$Lipkus - 5]

with(peters08[peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Positive frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
#abline(3.5,.625)

with(peters08[!peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Negative frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
#abline(3.5,.625)


###############################################
## Data creation
set.seed(1090)
surerisk <- list(y = peters08[,8:17], N = 108, J = 10)
sureriskreg <- list(y = peters08[,8:17], Frame = peters08$posframe,
                           Num = peters08$Lipkus, N = 108, J = 10)
choicesurerisk <- list(y = peters08[,8:17], N = 108, J = 10, K = 5,
                       X = peters08[,3:7],
                       Frame = peters08$posframe, 
                       Num = peters08$Lipkus,
                       Numc=mean(peters08$Lipkus))


###############################################
## Two-factor model estimation
if(file.exists("mod2res.rda")){
  load("mod2res.rda")
} else {
  mod2 <- jags.model("2f.jag", data = surerisk, n.chains = 3)
  update(mod2, 5000)
  mod2res <- coda.samples(mod2, c("lambda1", "lambda2","phi"), 
             n.iter=25000)
  save(mod2res, file="mod2res.rda")
}


###############################################
## Bayes factor for one- versus two-factor model
if(file.exists("logBF2f1f.rda")){
  load("logBF2f1f.rda")
} else {
  logBF2f1f <- logintlik(model = "2f.jag", 
                         paraest = c("lambda1", "lambda2", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_cov1"),
                         data = surerisk,
                         priorphi = "dwish") -
               logintlik(model = "1f.jag", 
                         paraest = c("lambda1", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_sig_f1"),
                         data = surerisk,
                         priorphi = "gamma")
  save(logBF2f1f, file="logBF2f1f.rda")
}  


###############################################
## Bayes factor for two-factor exploratory vs
## two-factor confirmatory model (not used in
## the paper)
if(file.exists("logBF2fc.rda")){
  load("logBF2fc.rda")
} else {
  logBF2fc <- logintlik(model = "2f.jag", 
                        paraest = c("lambda1", "lambda2", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_cov1"),
                        data = surerisk,
                        priorphi = "dwish") -
              logintlik(model = "2fre.jag", 
                        paraest = c("lambda1", "lambda2", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_cov1"),
                        data = surerisk,
                        prior = "dwish")
  save(logBF2fc, file="logBF2fc.rda")
}  


###############################################
## Table 1
loadings <- summary(mod2res)$statistics[1:20,1]
loadmat <- t(matrix(round(loadings, 2), 10, 2))
loadtab <- apply(loadmat, 1, function(x) paste(x, collapse=" & "))
loadtab <- paste(c("Riskless & ", "Risky & "), 
                 loadtab, "\\\\", sep="")
writeLines(loadtab)


## Frequentist measurement invariance tests
## (shhh, don't tell anybody this is here)
## ## Check what lavaan has to say
## library(lavaan)
## m1 <- 'f1 =~ Sure1.1 + Sure1.2 + Sure1.3 + Sure1.4 + Sure1.5
##        f2 =~ Risky1.1 + Risky1.2 + Risky1.3 + Risky1.4 + Risky1.5'
## m1fit <- cfa(m1, peters08, group="posframe")
## m4fit <- cfa(m1, peters08, group="posframe",
##              group.equal=c("loadings","intercepts", "means"))
## library(semTools)
## measurementInvariance(m1, peters08, group="posframe")


###############################################
## Measurement invariance Bayes factor (equal
## loadings or not)
if(file.exists("logBFmi1.rda")){
  load("logBFmi1.rda")
} else {
  ## loading are different wrt frame
 logBFmi1 <- logintlik(model = "mi_loading.jag", 
                       paraest = c("lambda1","lambda2",
                                   "delta1", "delta2", 
                                   "int1", "int2", "Inv_sig_e1", 
                                   "Inv_sig_e2", "Inv_cov1", 
                                   "Inv_cov2"),
                       data = c(surerisk, list(t=1)),
                       priorphi = "dwish") -
             logintlik(model = "mi_loading.jag", 
                       paraest = c("lambda1", "lambda2",
                                   "int1", "int2", "Inv_sig_e1", 
                                   "Inv_sig_e2", "Inv_cov1", 
                                   "Inv_cov2"),
                       data = c(surerisk, list(t=0)),
                       prior = "dwish")
  save(logBFmi1, file="logBFmi1.rda")
}  


###############################################
## Measurement invariance Bayes factor (equal
## intercepts or not (equal loadings assumed)
if(file.exists("logBFmi2.rda")){
  load("logBFmi2.rda")
} else {
  logBFmi2 <- logintlik(model = "mi_loadingintercept.jag", 
                        paraest = c("lambda1",
                                     "lambda2",
                                     "int1", "delta1", "delta2", 
                                     "delta3",
                                     "Inv_sig_e1", 
                                     "Inv_sig_e2", "Inv_cov1", 
                                     "Inv_cov2"),
                        data = c(surerisk, list(t=1)),
                        priorphi = "dwish") -
              logintlik(model = "mi_loadingintercept.jag", 
                        paraest = c("lambda1", "lambda2",
                                    "int1", "Inv_sig_e1", 
                                    "Inv_sig_e2", "Inv_cov1", 
                                    "Inv_cov2"),
                        data = c(surerisk, list(t=0)),
                        prior = "dwish")
  save(logBFmi2, file="logBFmi2.rda")
}


###############################################
## Two-factor model with covariates (no choice
## data included).  Not used in paper.
if(file.exists("summ2regres.rda")){
  load("summ2regres.rda")
} else {
    m2reg <- jags.model("2freg.jag", data = sureriskreg, 
                        n.chains = 3)
    update(m2reg,5000)
    m2regres <- coda.samples(m2reg, c("lambda1[1:5]",
                                      "lambda2[6:10]", 
                                      "beta0", "beta1", 
                                      "beta2", "beta3", "phi"),
                             n.iter=20000)
    summ2regres <- summary(m2regres)
    save(summ2regres, file="summ2regres.rda")
}
#summ2regres


###############################################
## Bayes factor for numeracy-by-frame interaction
## on riskless attraction
if(file.exists("logBFint.rda")){
  load("logBFint.rda")
} else {
  logBFint <- logintlik(model = "linkinteraction.jag", 
                        paraest = c("beta0", "lambda1", 
                                    "lambda2", 
                                    "beta1", "beta2", "beta3", 
                                    "Inv_sig_e1", "phi"),
                        data = c(sureriskreg, list(t=1)),
                        priorphi = "dunif") -
              logintlik(model = "linkinteraction.jag", 
                        paraest = c("beta0", "lambda1", 
                                    "lambda2", 
                                    "beta1", "beta2", 
                                    "Inv_sig_e1", "phi"),
                        data = c(sureriskreg, list(t=0)),
                        priorphi = "dunif")
  save(logBFint, file="logBFint.rda") 
}


###############################################
## Bayes factor for frame effect on riskless
## attraction
if(file.exists("logBFsure1.rda")){
  load("logBFsure1.rda")
} else {
  logBFsure1 <- logintlik(model = "linkFramesure1.jag", 
                         paraest = c("beta0", "lambda1", 
                                     "lambda2", "beta1", "beta2",
                                     "Inv_sig_e1", "phi"),
                         data = c(sureriskreg, list(t=1)),
                         priorphi = "dunif") -
                logintlik(model = "linkFramesure1.jag", 
                         paraest = c("beta0", "lambda1", 
                                     "lambda2", "beta2",
                                     "Inv_sig_e1", "phi"),
                         data = c(sureriskreg, list(t=0)),
                         priorphi = "dunif")
  save(logBFsure1, file="logBFsure1.rda")
}


###############################################
## Bayes factor for frame effect on risky
## attraction
if(file.exists("logBFrisky1.rda")){
  load("logBFrisky1.rda")
} else {
  logBFrisky1 <- logintlik(model = "linkFramerisky1.jag", 
                           paraest = c("beta0", "lambda1", 
                                       "lambda2", "beta1", "beta2",
                                       "Inv_sig_e1", "phi"),
                           data = c(sureriskreg, list(t=1)),
                           priorphi = "dunif") -
                 logintlik(model = "linkFramerisky1.jag", 
                           paraest = c("beta0", "lambda1", 
                                       "lambda2", "beta2",
                                       "Inv_sig_e1", "phi"),
                           data = c(sureriskreg, list(t=0)),
                           priorphi = "dunif")
  save(logBFrisky1, file="logBFrisky1.rda")
}


###############################################
## The main model, including two factors for
## attraction, choice data, numeracy, and so on
if(file.exists("summ2predres.rda")){
  load("summ2predres.rda")
} else {
    m2pred <- jags.model("2fpred.jag", data = choicesurerisk, 
                         n.chains = 3)
    update(m2pred,5000)
    m2predres <- coda.samples(m2pred, c("lambda1[1:5]", 
                                        "lambda2[6:10]", 
                                        "g1", "g2","g3","g4",
                                        "b0","b1","b2","b3",
                                        "beta0", "beta1", "beta2", 
                                        "beta3", "phi"), 
                              n.iter=20000)
    summ2predres <- summary(m2predres)
    save(summ2predres, file="summ2predres.rda")
}


###############################################
## Bayes factor for riskless attraction-by-numeracy
## interaction on choice
if(file.exists("logBFg4.rda")){
  load("logBFg4.rda")
} else {
  logBFg4 <- logintlik(model = "linkinteractionchoicesure.jag", 
                       paraest = c("beta0", "lambda1", 
                                   "lambda2", 
                                   "b0", "b1", "b2", "g1", "g2", 
                                   "g3","g4", "Inv_sig_ee",  
                                   "Inv_sig_e1", "phi"),
                       data = c(choicesurerisk, list(t=1)),
                       priorphi = "dunif") -
              logintlik(model = "linkinteractionchoicesure.jag", 
                        paraest = c("beta0", "lambda1", 
                                    "lambda2", 
                                    "b0", "b1", "b2", "g1", "g2", 
                                    "g3", "Inv_sig_ee",  
                                    "Inv_sig_e1", "phi"),
                        data = c(choicesurerisk, list(t=0)),
                        priorphi = "dunif")
  save(logBFg4, file="logBFg4.rda")
}


###############################################
## Bayes factor for risky attraction-by-numeracy
## interaction on choice
if(file.exists("logBFg3.rda")){
  load("logBFg3.rda")
} else {
  logBFg3 <- logintlik(model = "linkinteractionchoicerisky.jag", 
                       paraest = c("beta0", "lambda1", 
                                   "lambda2", 
                                   "b0", "b1", "b2", "g1", "g2", 
                                   "g3","g4", "Inv_sig_ee",  
                                   "Inv_sig_e1", "phi"),
                       data = c(choicesurerisk, list(t=1)),
                       priorphi = "dunif") -
              logintlik(model = "linkinteractionchoicerisky.jag", 
                        paraest = c("beta0", "lambda1", 
                                    "lambda2", 
                                    "b0", "b1", "b2", "g1", "g2", 
                                    "g4", "Inv_sig_ee",  
                                    "Inv_sig_e1", "phi"),
                        data = c(choicesurerisk, list(t=0)),
                        priorphi = "dunif")
  save(logBFg3, file="logBFg3.rda")
}


###############################################
## Bayes factor for the effect of frame on choice
if(file.exists("logBFfrchoice1.rda")){
  load("logBFfrchoice1.rda")
} else {
  logBFfrchoice1 <- logintlik(model = "linkFramechoice1.jag", 
                              paraest = c("beta0", "lambda1", 
                                          "lambda2", 
                                          "b0", "b1", "b2", 
                                          "g1", "g2", 
                                          "g3","g4",  
                                          "Inv_sig_ee",  
                                          "Inv_sig_e1", "phi"),
                           data = c(choicesurerisk, list(t=1)),
                                    priorphi = "dunif") -
                    logintlik(model = "linkFramechoice1.jag", 
                              paraest = c("beta0", "lambda1", 
                                          "lambda2", 
                                          "b0", "b2", "g1", "g2", 
                                          "g3", "g4", 
                                          "Inv_sig_ee",  
                                          "Inv_sig_e1", "phi"),
                              data = c(choicesurerisk, list(t=0)),
                                       priorphi = "dunif")
  save(logBFfrchoice1, file="logBFfrchoice1.rda")
}


###############################################
## Bayes factor for the frame-by-numeracy
## interaction on choice.
if(file.exists("logBFb3choice.rda")){
  load("logBFb3choice.rda")
} else {
  logBFb3choice <- logintlik(model = "linkb3choice.jag", 
                             paraest = c("beta0", "lambda1", 
                                         "lambda2", 
                                         "b0", "b1","b2", "b3",
                                         "g1", "g2", 
                                         "g3","g4",  "Inv_sig_ee",  
                                         "Inv_sig_e1", "phi"),
                             data = c(choicesurerisk, list(t=1)),
                             priorphi = "dunif") -
                   logintlik(model = "linkb3choice.jag", 
                             paraest = c("beta0", "lambda1", 
                                         "lambda2", 
                                         "b0", "b1", "b2", "g1", 
                                         "g2", 
                                         "g3", "g4", "Inv_sig_ee",  
                                         "Inv_sig_e1", "phi"),
                             data = c(choicesurerisk, list(t=0)),
                             priorphi = "dunif")
  save(logBFb3choice, file="logBFb3choice.rda")
}


###############################################
## Collect posterior intervals for Table 2;
## these are displayed with the above Bayes
## factors and F statistics from the original paper.
pname <- rownames(summ2predres$statistics)
rnums <- sapply(c("beta1", "b1", "beta3", "b3", "g3", "g4"), grep, 
                pname)
ints <- matrix(NA, 7, 2)
ints[1,] <- summ2predres$quantiles[rnums[[1]][1],c(1,5)]
ints[2,] <- summ2predres$quantiles[rnums[[1]][2],c(1,5)]
ints[3,] <- summ2predres$quantiles[rnums[[2]],c(1,5)]
ints[4,] <- summ2predres$quantiles[rnums[[3]],c(1,5)]
ints[5,] <- summ2predres$quantiles[rnums[[4]],c(1,5)]
ints[6,] <- summ2predres$quantiles[rnums[[5]],c(1,5)]
ints[7,] <- summ2predres$quantiles[rnums[[6]],c(1,5)]
ints <- apply(round(ints,2), 1, function(x) paste("(",paste(x, collapse=","),")",sep=""))


