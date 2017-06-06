
###############################################
## Preliminary packages, data, files
library("knitr")
## global knitr options
opts_chunk$set(echo=FALSE, fig.path='figure/vuong-',fig.align='center',fig.show='hold',size='footnotesize', cache.path="cache/", warning=FALSE, message = FALSE)

## packages
library("rjags", quietly=TRUE)
library("mvtnorm")
library("MCMCpack")
library("Matrix")
library("reshape2")
library("ggplot2")
library("polspline") ## for S-D bayes

## auxiliary code
source("bfcalc.R")
source("jagsmodels.R")

## real-world data
load("peters08.rda")


###############################################
## Figure 2

## correlations between observed vars
pfcor <- round(cor(peters08[peters08$posframe,c(8:17,3:7)]),2)
nfcor <- round(cor(peters08[!peters08$posframe,c(8:17,3:7)]),2)
pfcor[upper.tri(pfcor)] <- nfcor[upper.tri(nfcor)]
rownames(pfcor) <- colnames(pfcor) <- c(paste("AR",1:10,sep=""),paste("C",1:5,sep=""))
meltpfcor <- melt(pfcor[1:15,15:1])
ggplot(meltpfcor, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + geom_text(aes(Var2, Var1, label=value), size=2.5) + 
    scale_fill_gradient(name=expression("Pearson r"), low = "#fdf6e3", high = "steelblue",
    breaks=seq(-.3, 1, by = 0.2), limits = c(-.3, 1)) + xlab("") + ylab("") +
    theme(axis.ticks=element_blank(), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7))


###############################################
## Some other data visualizations that didn't
## make it into the paper.

## par(mfrow=c(1,2))
## shades <- paste("grey", c(100, 95, 85, 70, 60, 30))
## ptshading <- shades[peters08$Lipkus - 5]
## 
## with(peters08[peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Positive frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
## #abline(3.5,.625)
## 
## with(peters08[!peters08$posframe,], plot(jitter(Risky1.1 - Sure1.1, factor=1.5), jitter(sweden), pch=21, bg=ptshading, xlab="Risky rating - Riskless rating", ylab="Choice rating", main="Negative frame", xlim=c(-4.2,4.2), ylim=c(.8, 7.2)))
## #abline(3.5,.625)


## long1 <- melt(peters08, id=c(1,28:37), measure=8:17)
## long1$scenario <- rep(rep(1:5, each=108), 2)
## long1$riskless <- rep(c("Riskless", "Risky"), each=108*5)
## long3 <- melt(peters08, id=c(1,28:37), measure=3:7)
## petlong <- cbind(long1, c(long3[,13], long3[,13]))
## names(petlong)[c(13,16)] <- c("attraction", "choice")
## petlong$frame <- rep("Negative", nrow(petlong))
## petlong$frame[petlong$posframe] <- "Positive"
## 
## ggplot(petlong, aes(x=attraction, y=choice, colour=Lipkus)) + facet_wrap(~ scenario + frame + riskless, nrow=5) + geom_point(position=position_jitter(w=.2, h=.2)) + xlab("Attraction Rating") + ylab("Choice Rating")
## 
## petlong2 <- petlong[grep("Sure", petlong$variable),]
## petlong2$riskyatt <- petlong[grep("Risky", petlong$variable),]$attraction
## ggplot(petlong2, aes(x=attraction, y=riskyatt)) + facet_wrap(~ scenario + frame, nrow=5) + geom_point(position=position_jitter(w=.2, h=.2)) + xlab("Riskless Attraction") + ylab("Risky Attraction")


###############################################
## Data creation

## col 8:12 represents attractiveness rating of sure thing.
## col 13:17 represents attractiveness rating of risky thing.
## N: the number of observations. 
## J: 2*the number of scenarios(5)
## get all data set used later. 
set.seed(1090)
surerisk <- list(y = peters08[,8:17], N = 108, J = 10, halfJ = 10/2)
sureriskreg <- list(y = peters08[,8:17], Frame = peters08$posframe,
                           Num = peters08$Lipkus, N = 108, J = 10,
                    halfJ = 10/2)
choicesurerisk <- list(y = peters08[,8:17], N = 108, J = 10, K = 5,
                       X = peters08[,3:7],
                       Frame = peters08$posframe, 
                       Num = peters08$Lipkus,
                       Numc= mean(peters08$Lipkus),
                       halfJ = 10/2)


###############################################
## Two-factor model estimation
if(file.exists("mod2res.rda")){
  load("mod2res.rda")
} else {
  mod2 <- jags.model("2f.jag", data = surerisk, 
                     n.chains = 3)
  update(mod2, 5000)
  mod2res <- coda.samples(mod2, c("lambda1", "lambda2","phi"), 
             n.iter=5000)
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
                         priorphi = "dwish", n.iter = 5000) -
               logintlik(model = "1f.jag", 
                         paraest = c("lambda1", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_sig_f1"),
                         data = surerisk,
                         priorphi = "gamma", n.iter = 5000)
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
                        priorphi = "dwish", n.iter = 5000) -
              logintlik(model = "2fre.jag", 
                        paraest = c("lambda1", "lambda2", 
                                     "int1", "Inv_sig_e1",
                                     "Inv_cov1"),
                        data = surerisk,
                        prior = "dwish", n.iter = 5000)
  save(logBF2fc, file="logBF2fc.rda")
}  


###############################################
## Table 1
loadings <- summary(mod2res)$statistics[1:20,1]
loadmat <- t(matrix(formatC(loadings, 2, format='f'), 10, 2))
loadmat[1,6] <- "0"
loadmat[2,1] <- "0"
loadtab <- apply(loadmat, 1, function(x) paste(x, collapse=" & "))
loadtab <- paste(c("Riskless & ", "Risky & "), 
                 loadtab, "\\\\", sep="")
writeLines(loadtab)


###############################################
## Functions for specifying prior precisions
para <- c("beta0", "lambda1","lambda2", "lambda3", 
          "b0", "beta1[1]","beta1[2]", "beta2[1]","beta2[2]", 
          "beta3[1]", "beta3[2]","beta4","beta5",
          "beta6", "nu1", "nu2","nu3", "nu4",
          "Inv_sig_ee", "Inv_sig_e1", "phi")
p <- length(para)
tlist <- function(tbeta1_1=1, tbeta1_2=1, tbeta2_1=1, tbeta2_2=1,  
                  tbeta3_1=1, tbeta3_2=1, tbeta4=1, tbeta5=1, 
                  tbeta6=1, tnu1=1, tnu2=1, tnu3=1, tnu4=1){
           list(tbeta1_1 = tbeta1_1, 
                tbeta1_2 = tbeta1_2,
                tbeta2_1 = tbeta2_1,
                tbeta2_2 = tbeta2_2,
                tbeta3_1 = tbeta3_1,
                tbeta3_2 = tbeta3_2,
                tbeta4 = tbeta4, 
                tbeta5 = tbeta5, 
                tbeta6 = tbeta6, 
                tnu1 = tnu1, tnu2 = tnu2,
                tnu3 = tnu3, tnu4 = tnu4)
         }
prlist <- function(plambda1=.1, plambda2=.1, plambda3=.1, 
                  pbeta1=.1, pbeta2=.1, pbeta3=.1,  
                  pbeta4=.1, pbeta5=.1, 
                  pbeta6=.1, pnu1=.1, pnu2=.1, pnu3=.1, 
                  pnu4=.1){
           list(plambda1 = plambda1, 
                plambda2 = plambda2,
                plambda3 = plambda3,
                pbeta1 = pbeta1, 
                pbeta2 = pbeta2,
                pbeta3 = pbeta3,
                pbeta4 = pbeta4, 
                pbeta5 = pbeta5, 
                pbeta6 = pbeta6, 
                pnu1 = pnu1, pnu2 = pnu2,
                pnu3 = pnu3, pnu4 = pnu4)
       }


###############################################
## Fitting the SEM from Figure 1
if(file.exists("semresp.rda")){
  load("semresp.rda")
} else {
  if(!file.exists("semresp.rda")){
    plist <- prlist()
    sem <- jags.model("2fsemt.jag", c(choicesurerisk, tlist(), 
                                      plist), n.chains = 3)
    update(sem, 5000)
    semresp <- coda.samples(sem, para, n.iter=5000)    
    save(semresp, file="semresp.rda")
  }
}


###############################################
## Laplace log-Bayes factors from the SEM
if(file.exists("logBFsemp.rda")){
  load("logBFsemp.rda")
} else {
  ## full model integrated log likelihood
  band <- 1
  npara <- 13
  logintsem <- matrix (NA, band, npara)
  full <- rep(NA, nrow(logintsem))
  plist <- prlist()
  for (i in 1:nrow(logintsem)){
    full[i] <- logintlik(model = "2fsemt.jag", 
                         paraest = para,
                         data = c(choicesurerisk, tlist(), 
                         plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta1_1(riskless--frame)
    logintsem[i,1] <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta1[1]",p))],
                       data = c(choicesurerisk, tlist(tbeta1_1=0), 
                                plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta1_2(risky--frame)
    logintsem[i,2]  <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta1[2]",p))],
                       data = c(choicesurerisk, tlist(tbeta1_2=0), 
                                plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta2_1(riskless--Num)
    logintsem[i,3]  <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta2[1]",p))],
                       data = c(choicesurerisk, tlist(tbeta2_1=0), 
                                plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta2_2(risky--Num)
    logintsem[i,4]  <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta2[2]",p))],
                       data = c(choicesurerisk, tlist(tbeta2_2=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta3_1(riskless--Num*frame)
    logintsem[i,5]  <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta3[1]",p))],
                       data = c(choicesurerisk, tlist(tbeta3_1=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta3_2(risky--Num*frame)
    logintsem[i,6]  <- logintlik(model = "2fsemt.jag", 
                   paraest = para[-which(para==rep("beta3[2]",p))],
                       data = c(choicesurerisk, tlist(tbeta3_2=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta4(choice--frame)
    logintsem[i,7]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("beta4",p))],
                         data = c(choicesurerisk, tlist(tbeta4=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta5(choice--Num)
    logintsem[i,8]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("beta5",p))],
                         data = c(choicesurerisk, tlist(tbeta5=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## beta6(choice--Num*frame interaction)
    logintsem[i,9]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("beta6",p))],
                         data = c(choicesurerisk, tlist(tbeta6=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## nu1(choice--riskless)
    logintsem[i,10]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("nu1",p))],
                         data = c(choicesurerisk, tlist(tnu1=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## nu2(choice--risky)
    logintsem[i,11]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("nu2",p))],
                         data = c(choicesurerisk, tlist(tnu2=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## nu3(choice--riskless*Num interaction)
    logintsem[i,12]  <- logintlik(model = "2fsemt.jag", 
                      paraest = para[-which(para==rep("nu3",p))],
                         data = c(choicesurerisk, tlist(tnu3=0), plist), 
                         priorphi = "dunif", n.iter=5000)
  ## nu4(choice--risky*Num interaction)
    logintsem[i,13]  <- logintlik(model = "2fsemt.jag", 
                        paraest = para[-which(para==rep("nu4",p))],
                        data = c(choicesurerisk, tlist(tnu4=0), plist), 
                        priorphi = "dunif", n.iter=5000)
    }
    ## two population sd. (based on t test). 
    sdpool <- sqrt((rep(var(full), ncol(logintsem))+
                   apply(logintsem, 2, sd)^2)/nrow(logintsem))
    ## combine mean difference, pooled sd. 95% lower bound, 
    ## 95% upper bound
    logBFsem <- cbind(mean(full)-(colMeans(logintsem)), sdpool, 
                ((mean(full)-colMeans(logintsem)) - qt(.975, 
                ((apply(logintsem, 2, sd)^2+var(full))^2*
                     nrow(logintsem)-1)/(
               apply(logintsem, 2, sd)^4 + var(full)^2))*sdpool), 
                ((mean(full)-colMeans(logintsem)) + qt(.975, 
                 ((apply(logintsem, 2, sd)^2+var(full))^2*
                     nrow(logintsem)-1)/(
               apply(logintsem, 2, sd)^4 + var(full)^2))*sdpool))
  
    rownames(logBFsem) <- c("logBFbeta1_1", "logBFbeta1_2", 
                       "logBFbeta2_1", "logBFbeta2_2", 
                       "logBFbeta3_1", "logBFbeta3_2",          
                       "logBFbeta4", "logBFbeta5", "logBFbeta6",
                       "logBFnu1","logBFnu2", "logBFnu3", 
                       "logBFnu4")
    colnames(logBFsem) <- c("mean","sdpool","95low","95up")
    if(!file.exists("logBFsemp.rda")){
      logBFsemp <- logBFsem
      save(logBFsemp, file="logBFsemp.rda")}
}


###############################################
## Savage-Dickey log-Bayes factors for our
## comparison (not included in paper)
## ## non-informative results
## # res <- rbind(semres[[1]], semres[[2]], semres[[3]])
## ## informative results
##  res <- rbind(semresp[[1]], semresp[[2]], semresp[[3]])
## 
## ## beta1_1(riskless--frame).
## beta1_1 <- as.data.frame(res)$"beta1[1]"
## postdens <- logspline(beta1_1)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[4])))
## logBF01beta1_1 <- log(postest/priorest)
## 
## ## beta1_2(risky--frame).
## beta1_2 <- as.data.frame(res)$"beta1[2]"
## postdens <- logspline(beta1_2)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[4])))
## logBF01beta1_2 <- log(postest/priorest)
## 
## ## beta2_1(riskless--Num)
## beta2_1 <- as.data.frame(res)$"beta2[1]"
## postdens <- logspline(beta2_1)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[5])))
## logBF01beta2_1 <- log(postest/priorest)
## 
## ## beta2_2(risky--Num)
## beta2_2 <- as.data.frame(res)$"beta2[2]"
## postdens <- logspline(beta2_2)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[5])))
## logBF01beta2_2 <- log(postest/priorest)
## 
## ## beta3_1(riskless--frame*Num interaction).
## beta3_1 <- as.data.frame(res)$"beta3[1]"
## postdens <- logspline(beta3_1)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[6])))
## logBF01beta3_1 <- log(postest/priorest)
## 
## ## beta3_2(risky--frame*Num interaction)
## beta3_2 <- as.data.frame(res)$"beta3[2]"
## postdens <- logspline(beta3_2)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[6])))
## logBF01beta3_2 <- log(postest/priorest)
## 
## ## beta4(choice--frame).
## beta4 <- as.data.frame(res)$beta4
## postdens <- logspline(beta4)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[7])))
## logBF01beta4 <- log(postest/priorest)
## 
## ## beta5(choice--Num)
## beta5 <- as.data.frame(res)$beta5
## postdens <- logspline(beta5)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[8])))
## logBF01beta5 <- log(postest/priorest)
## 
## ## beta6(choice--frame*Num interaction).
## beta6 <- as.data.frame(res)$beta6
## postdens <- logspline(beta6)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[9])))
## logBF01beta6 <- log(postest/priorest)
## 
## ## nu1(choice--riskless).
## nu1 <- as.data.frame(res)$nu1
## postdens <- logspline(nu1)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[10])))
## logBF01nu1 <- log(postest/priorest)
## 
## ## nu2(choice--risky).
## nu2 <- as.data.frame(res)$nu2
## postdens <- logspline(nu2)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(prlist()[11])))
## # priorest <- dunif(0, -1, 1)
## logBF01nu2 <- log(postest/priorest)
## 
## ## nu3(choice--riskless*Num interaction).
## nu3 <- as.data.frame(res)$nu3
## postdens <- logspline(nu3)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[12])))
## logBF01nu3 <- log(postest/priorest)
## 
## ## Validate with SD nu4(choice--risky*Num interaction).
## nu4 <- as.data.frame(res)$nu4
## postdens <- logspline(nu4)
## postest <- dlogspline(0,postdens)
## priorest <- dnorm(0,0,sqrt(1/as.numeric(plist[13])))
## logBF01nu4 <- log(postest/priorest)
## 
## logBF01 <- c(logBF01beta1_1, logBF01beta1_2, logBF01beta2_1,
##              logBF01beta2_2, logBF01beta3_1, logBF01beta3_2,
##              logBF01beta4, logBF01beta5,
##              logBF01beta6, logBF01nu1, logBF01nu2,
##              logBF01nu3, logBF01nu4)


###############################################
## Collect posterior intervals for Table 2;
## these are displayed with the above log-Bayes
## factors and F statistics from the original paper.
pname <- rownames(summary(semresp)$statistics)
rnums <- sapply(c("beta1[1]", "beta1[2]", "beta2[1]", "beta2[2]", 
                  "beta3[1]", "beta3[2]", 
                  "beta4", "beta5", "beta6", "nu1","nu2", 
                  "nu3","nu4"), grep, 
                pname, fixed=TRUE)

intsp <- summary(semresp)$quantiles[rnums,c(1,5)]
intsp <- apply(round(intsp,2), 1, function(x) paste("(",paste(x, collapse=","),")",sep=""))


