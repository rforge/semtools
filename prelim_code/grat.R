################
## Prelims
library("foreign")
library("lavaan")
library("strucchange")
## I assume these are in your path:
source("efpFunctional-cat.R")
source("lavaan-S3.R")

info_full <- function(x, ...) solve(vcov(x) * nobs(x))

## The warning here is meaningless, I think:
dat <- read.spss("grat data combined_9.29.09.sav", to.data.frame=TRUE)
names(dat) <- tolower(names(dat))

################
## Data notes:
## Not sure why AO_4 has some above the maximum of 9; remove them
dat <- dat[dat$ao_4 <= 9 ,]

## Not sure why some people have scale values that are not whole numbers
apply(dat[,4:28],2,function(x){sum((x - floor(x)) > 0)})
## Round these for now (but Fan is checking on this)
dat[,4:28] <- round(dat[,4:28])

## GRAT Item 2 excluded, which is losd_1 in the data.
## gq6_6 is excluded; see bottom p 313 to top p 314

## Approximate replication of one of their models; could
## continue this for all of Table 1.  The test statistics are
## generally off, but the loadings are pretty close.  Note that,
## when they do the measurement invariance tests following Table 1, they 
## set "std.lv=FALSE".
m1 <- cfa('f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
          data=dat[dat$agegroup=="12-13 yr olds",],
          std.lv=TRUE, std.ov=TRUE)
summary(m1)

################
## Some models

## Fit a one-group model for GQ-6 and test for invariance
## in loadings:
m2 <- cfa('f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
          data=dat, meanstructure=TRUE)

gefp_2 <- gefp(m2, fit = NULL, scores=estfun.lavaan,
               order.by = as.numeric(dat$agegroup), 
               vcov = info_full, sandwich = FALSE, parm=1:4)
sctest(gefp_2, functional=catL2BB(gefp_2))  # p = .01
sctest(gefp_2, functional=ordL2BB(gefp_2))  # p = .1

## Fit a one-group model for GAC and test for invariance
## in loadings:
m3 <- cfa('f1 =~ gac_1 + gac_2 + gac_3',
          data=dat, meanstructure=TRUE)

gefp_3 <- gefp(m3, fit = NULL, scores=estfun.lavaan,
               order.by = as.numeric(dat$agegroup), 
               vcov = info_full, sandwich = FALSE, parm=1:2)
sctest(gefp_3, functional=catL2BB(gefp_3))  # p = .67
sctest(gefp_3, functional=ordL2BB(gefp_3))  # p = .35


## Fit a one-group model for GRAT and test for invariance
## (GRAT Item 2 excluded, which is losd_1 in the data)
m4 <- cfa('f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
           f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
           f3 =~ ao_1 + ao_2 + ao_3 + ao_4
           f1 ~ 0*1
           f2 ~ 0*1
           f3 ~ 0*1', data=dat,
           meanstructure=TRUE)

gefp_m4 <- gefp(m4, fit = NULL, scores=estfun.lavaan,
                order.by = as.numeric(dat$agegroup), 
                vcov = info_full, sandwich = FALSE, parm=1:4)
## The catL2BB p-value is sometimes non-significant, depending on
## how one handles the errors in the data.  Also might get something
## by changing the identification constraint (choose a different
## factor loading to fix at one).
sctest(gefp_m4, functional=catL2BB(gefp_m4))  # p = .047
ordfun <- ordL2BB(gefp_m4)
sctest(gefp_m4, functional=ordfun)  # p = .001
## The plot has five categories, but there are actually six defined by
## dat$agegroup.  I haven't looked at the code that defines the ordinal
## test statistic, so maybe this is good behavior.
plot(gefp_m4, functional=ordfun)
abline(h=ordfun$computeCritval(alpha=.05, nproc=4), col=4)

## Could also test other sets of parameters, define age groups differently
## (could use each year on its own, but time to compute critical value
## increases greatly), or test within
## multiple group models where some parameters vary from group to
## group and some parameters are fixed across groups.  However, I think
## one of the cool things about the ordinal test is that it does not
## immediately correspond to a LRT with multiple groups (or maybe it
## corresponds to a LRT with multiple groups, where there are order
## constraints in the parameters?).
