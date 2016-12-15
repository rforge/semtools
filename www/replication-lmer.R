###################################################
### Preliminaries
###################################################

## packages
library("lavaan")
library("lme4")
library("reshape2")
library("xtable")
library("sandwich")
library("data.table")

## code
source("estfun.lmerMod.R")
source("vcov.full.lmerMod.R")
source("bread.lmerMod.R")


###################################################
### Model fit via lme4
###################################################

library("lme4")
lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)


###################################################
### Equivalent model fit via lavaan
###################################################

## data format manipulation
testwide <- reshape2::dcast(sleepstudy, Subject ~ Days, value.var = "Reaction")
names(testwide)[2:11] <- paste("d", 1:10, sep="")

## model fit      
latent <- '
      i =~ 1*d1 + 1*d2 + 1*d3 + 1*d4 + 1*d5
           + 1*d6 + 1*d7 + 1*d8 + 1*d9 + 1*d10
      
      s = ~ 0*d1 + 1*d2 + 2*d3 + 3*d4 + 4*d5
            + 5*d6 + 6*d7 + 7*d8 + 8*d9 + 9*d10

      d1 ~~ evar*d1
      d2 ~~ evar*d2
      d3 ~~ evar*d3
      d4 ~~ evar*d4
      d5 ~~ evar*d5
      d6 ~~ evar*d6
      d7 ~~ evar*d7
      d8 ~~ evar*d8
      d9 ~~ evar*d9
      d10 ~~ evar*d10

      ## reparameterize as sd
      sdevar := sqrt(evar)
      i ~~ ivar*i
      isd := sqrt(ivar)
      '
lavaanfit <- growth(latent, data = testwide, estimator="ML")


###################################################
### Parameter estimates comparison
###################################################

rbind(round(coef(lavaanfit)[c(14:15,11:13,1)], 2),
      round(c(fixef(lme4fit),
              as.numeric(VarCorr(lme4fit)$Subject)[c(1,4,3)],
              getME(lme4fit, "sigma")^2), 2))


###################################################
### Gradients
###################################################

## casewise
score1 <- estfun.lmerMod(lme4fit, level = 1)
gradients1 <- colSums(score1)
gradients1

## clusterwise
score2 <- estfun.lmerMod(lme4fit, level= 2)
gradients2 <- colSums(score2)
gradients2


###################################################
### Figure 1
###################################################

score2_lavaan <- estfun.lavaan(lavaanfit)
score2_lavaan <- as.data.table(score2_lavaan)
score2_lavaan <- setcolorder(score2_lavaan, c(5, 6, 2, 4, 3, 1))

plot(as.matrix(score2_lavaan), as.matrix(score2), 
     xlab="estfun.lavaan", 
     ylab = "estfun.lmerMod", pch=20)


###################################################
### Table 1
###################################################

test <- vcov(lavaanfit, remove.duplicated=TRUE)
## reorder
test <- as.data.table(test)
test <- setcolorder(test, c(5, 6, 2, 4, 3, 1))
orderlme4 <- c(6, 3, 5, 4, 1, 2)
test <-  test[order(as.numeric(orderlme4)),]
test <- round(test, 2)
## create tab
tab <- matrix(NA, 36, 5)
tab[,3:4] <- format(round(cbind(as.vector(vcov.full.lmerMod(lme4fit)), 
                         as.vector(as.matrix(test))),2),nsmall=2)
tab[,1] <- rep(colnames(bread(lme4fit)), 6)
tab[,2] <- rep(colnames(bread(lme4fit)), each=6)
tab[,5] <- format(round(abs(as.numeric(tab[,3])-as.numeric(tab[,4])),2),
                  nsmall=2)
colnames(tab) <- c("Column Name", "Row Name", "Bread", "Lavaan", "Abs(diff)")

## print latex tab by using xtable()
print(xtable(tab, caption="Comparison between \\code{vcov.full.lmerMod()} 
                           output and \\pkg{lavaan} \\code{vcov()} output for 
                           the \\code{SleepStudy} model. The first two columns 
                           describe the specific matrix entry being compared, 
                           the third and fourth columns show the estimates, 
                           and the fifth column shows the absolute difference.", 
                           label="tab:bread"), 
      include.rownames=FALSE)


###################################################
### Table 2
###################################################

sand <- sandwich(lme4fit, bread. = bread.lmerMod(lme4fit), 
                 meat. = meat(lme4fit, level = 2))

## directly produced by lavaan
lavaanfitrob <- growth(latent, data = testwide, estimator = "MLR")
test2 <- vcov(lavaanfitrob, remove.duplicated=TRUE)
## reorder 
test2 <- as.data.table(test2)
test2 <- setcolorder(test2, c(5, 6, 2, 4, 3, 1))
orderlme4 <- c(6, 3, 5, 4, 1, 2)
test2 <-  test2[order(as.numeric(orderlme4)),]
test2 <- round(test2, 2)
## create tab
tab2 <- matrix(NA, 36, 5)
tab2[,3:4] <- format(round(cbind(as.vector(sand), 
                                 as.vector(as.matrix(test2))),2), nsmall=2)
tab2[,1] <- rep(colnames(sand), 6)
tab2[,2] <- rep(colnames(sand), each=6)
tab2[,5] <- format(round(abs(as.numeric(tab2[,3])-as.numeric(tab2[,4])),2),
                   nsmall=2)
colnames(tab2) <- c("Column Name", "Row Name", "Bread", "Lavaan", "Abs(diff)")

## print latex tab by using xtable()
print(xtable(tab2,caption="Comparison of the \\code{SleepStudy} sandwich 
                           estimator obtained from our \\code{lmerMod} code 
                           with the analogous estimator obtained from 
                           \\pkg{lavaan}. The first two columns describe the 
                           specific matrix entry being compared, the third and 
                           fourth columns show the estimates, and the fifth 
                           column shows the absolute difference.", 
             label="tab:sandwich"), include.rownames=FALSE)
