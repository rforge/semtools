## Code to replicate examples from Merkle, Furr, Rabe-Hesketh

## convenience functions:
source("mfrh_functions.R")

######################
## CFA Example
## requires blavaan 0.3-2 or higher; download at
## install.packages("blavaan", repos="http://faculty.missouri.edu/~merklee", type="source")
library("blavaan")
library("coda")
## data
library("psychotools")

data("StereotypeThreat", package="psychotools")
StereotypeThreat <- transform(StereotypeThreat, group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]

## hold results
niter <- 10
res <- vector("list", niter)
npar <- rep(NA, 9)
for(i in 1:niter){
  res[[i]] <- list(dicres = matrix(NA, 9, 8),
                   llres = matrix(NA, 9, 4),
                   pdres = matrix(NA, 9, 8))
}

## MGCFA
## convenience function for multi-group CFA on this data
mgcfa <- function(model, ...) bcfa(model, data = StereotypeThreat,
                                   group = "group", burnin=4000,
                                   sample=1e4, save.lvs=TRUE, jags.ic=TRUE,
                                   dp = dpriors(ipsi="dgamma(1,1)"),
                                   inits = "simple", ...)

for(i in 1:niter){
  ## Step 2: Fix loadings across groups
  f <- 'ability =~ abstract + verbal + numerical'
  m2 <- mgcfa(f, group.equal = "loadings")
    
  res <- fillres(m2, res, 1, i)

  ## Step 2a: Free numerical loading in group 4 (minority.threat)
  f <- 'ability =~ abstract + verbal + c(l1, l1, l1, l4) * numerical'
  m2a <- mgcfa(f, group.equal = "loadings")
  
  res <- fillres(m2a, res, 2, i)

  ## Step 3: Fix variances across groups
  m3 <- mgcfa(f, group.equal = c("loadings", "residuals"))
  
  res <- fillres(m3, res, 3, i)

  ## Step 3a: Free numerical variance in group 4
  f <- c(f, 'numerical ~~ c(e1, e1, e1, e4) * numerical')
  m3a <- mgcfa(f, group.equal = c("loadings", "residuals"))

  res <- fillres(m3a, res, 4, i)
    
  ## Step 4: Fix latent variances within conditions
  f <- c(f, 'ability ~~ c(vmaj, vmin, vmaj, vmin) * ability')
  m4 <- mgcfa(f, group.equal = c("loadings", "residuals"))

  res <- fillres(m4, res, 5, i)
    
  ## Step 5: Fix certain means, free others
  f <- c(f, 'numerical ~ c(na1, na1, na1, na4) * 1')
  m5 <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

  res <- fillres(m5, res, 6, i)
    
  ## Step 5a: Free ability mean in group majority.control
  f <- c(f, 'abstract ~ c(ar1, ar2, ar2, ar2) * 1')
  m5a <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

  res <- fillres(m5a, res, 7, i)
    
  ## Step 5b: Free also ability mean in group minority.control
  f <- c(f[1:4], 'abstract ~ c(ar1, ar2, ar3, ar3) * 1')
  m5b <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))
  
  res <- fillres(m5b, res, 8, i)
    
  ## Step 6: Different latent mean structure
  f <- c(f, 'ability ~  c(maj, min1, maj, min2) * 1 + c(0, NA, 0, NA) * 1')
  m6 <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

  res <- fillres(m6, res, 9, i)
}
save(res, file="semdicres3.rda")
  
npar[1] <- fitMeasures(m2, "npar")
npar[2] <- fitMeasures(m2a, "npar")
npar[3] <- fitMeasures(m3, "npar")
npar[4] <- fitMeasures(m3a, "npar")
npar[5] <- fitMeasures(m4, "npar")
npar[6] <- fitMeasures(m5, "npar")
npar[7] <- fitMeasures(m5a, "npar")
npar[8] <- fitMeasures(m5b, "npar")
npar[9] <- fitMeasures(m6, "npar")
save(npar, file="semnpar.rda")


#############
## Figure 1
#############
justdics <- lapply(res, function(x) x$dicres)
alldics <- do.call("rbind", justdics)
alldics <- cbind.data.frame(alldics, model=rep(1:9, 10))
names(alldics)[1:8] <- c("dic", "dicse", "jdic", "jdicse",
                         "conddic", "conddicse", "condjdic",
                         "condjdicse")
dic <- as.numeric(as.matrix(alldics[,c(1,3,5,7)]))
dicse <- as.numeric(as.matrix(alldics[,c(2,4,6,8)]))
mod <- as.factor(rep(c("2","2a","3","3a","4","5","5a","5b","6"), 4))
type <- rep(rep(c("Spiegelhalter","Plummer"), each=nrow(alldics)), 2)
marg <- rep(c("Marginal","Conditional"), each=2*nrow(alldics))
plotdic <- cbind.data.frame(dic, dicse, mod, type, marg)

fakedata <- data.frame(type=rep("Spiegelhalter",4), marg=rep(c("Marginal","Conditional"), each=2), mod=rep("2",4), dic=c(4307,4332,3959,4209))

my_theme <- theme_bw() +
  theme(text = element_text(family = "serif", size = 14),
        axis.text.x = element_text(size = 11),
        strip.background = element_rect(fill = NA, color = NA, size = .5),
        legend.key = element_rect(fill = NA, color = NA))

# NB to get JAGS error bars, use type=="JAGS" and uncomment last
# part
ggplot(subset(plotdic, type=="Spiegelhalter")) +
  aes(y = dic, x = as.factor(mod)) +
  geom_jitter(height = 0, width = .1, alpha = .25) +
  facet_wrap(~ marg, scale = "free_y") +
  geom_blank(data=fakedata) +
  labs(y = "DIC value", x = "Model") +
  my_theme +
  theme(panel.grid.minor.x = element_blank())# + geom_errorbar(aes(ymin=dic-dicse,ymax=dic+dicse), width=.1)


#############
## Figure 2
#############
justpd <- lapply(res, function(x) x$pdres)
allpd <- do.call("rbind", justpd)
allpd <- cbind.data.frame(allpd, model=rep(1:9, 10))
names(allpd)[1:8] <- c("dic", "dicse", "jdic", "jdicse",
                       "conddic", "conddicse", "condjdic",
                       "condjdicse")
pd <- as.numeric(as.matrix(allpd[,c(1,3,5,7)]))
mod <- as.factor(rep(c("2","2a","3","3a","4","5","5a","5b","6"), 4))
#mod <- rep(allpd$model, 4)
type <- rep(rep(c("Spiegelhalter","Plummer"), each=nrow(allpd)), 2)
marg <- rep(c("Marginal","Conditional"), each=2*nrow(allpd))
plotpd <- cbind.data.frame(pd, mod, type, marg, npar)

library(reshape2)
ppdnew <- melt(plotpd, id.vars = c('mod', 'type', 'marg'))
cnrows <- which(ppdnew$variable == "npar" & ppdnew$marg == "Conditional")
ppdnew$value[cnrows] <- ppdnew$value[cnrows] + 295 - rep(nlvpar, length(cnrows)/9)

ggplot(ppdnew, aes(x = as.factor(mod), y = value)) +
  geom_jitter(data = ppdnew[ppdnew$variable == 'pd',], height = 0, width = .1, alpha = .25) +
  geom_line(aes(x = as.numeric(mod)), data = ppdnew[ppdnew$variable == 'npar',]) + 
  facet_grid(marg ~ type, scale = "free_y") +
  labs(y = "pD value", x = "Model") +
  my_theme +
  theme(panel.grid.minor.x = element_blank())


#############
## Figure 3
#############
jagdic <- subset(plotdic, type=="Plummer")
jagdic$rep <- rep(rep(1:10, each=9), 2)
delreps <- which(jagdic$rep != 6)
jagdic$dicse[delreps] <- NA

ggplot(jagdic) +
  aes(y = dic, x = as.factor(mod)) +
  geom_jitter(height = 0, width = .1, alpha = .25) +
  facet_wrap(~ marg, scale = "free_y") +
  geom_blank(data=fakedata) +
  labs(y = "DIC value", x = "Model") +
  my_theme +
  theme(panel.grid.minor.x = element_blank()) + geom_errorbar(aes(ymin=dic-2*dicse,ymax=dic+2*dicse), width=.1)



######################
## IRT Example
##
library(rstan)
library(edstan)
library(loo)
library(reshape2)
library(doParallel)
options(mc.cores = 5)
options(loo.cores = 5)

# Assemble example dataset
dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = ~ 1 + male + anger)

# Fit model
fit <- stan("rasch_edstan_modified.stan", data = dl, iter = 500, chains = 5)

# Obtain marginal likelihoods
cl <- makeCluster(5)
registerDoParallel(cl)
ll_marg <- mll_parallel(fit, dl, f_marginal, "zeta", "sigma", 11)
stopCluster(cl)

# Obtain marginal information criteria
dic(ll_marg)
waic(ll_marg$ll)
loo(ll_marg$ll)
