###################
## Packages/code ##
###################

## packages
library("foreign")
library("lavaan")
library("strucchange")

## additional code
source("efpFunctional-cat.R")
source("../www/lavaan-S3.R")
info_full <- function(x, ...) solve(vcov(x) * nobs(x))


##########
## Data ##
##########

## with warnings (but probably harmless)
YouthGratitude <- read.spss("grat_data_combined_2009-09-29.sav",
  to.data.frame = TRUE)

## recode/unify variable names and factor levels
names(YouthGratitude) <- tolower(names(YouthGratitude))
levels(YouthGratitude$agegroup) <- sapply(strsplit(levels(YouthGratitude$agegroup), " "), head, 1)

## reorder columns to match order of usage in the paper
YouthGratitude <- YouthGratitude[, c(1:3, 20:28, 4:19)]

## data processing required due to scale values that are not
## integers or outside of range, hence make copy
yg <- YouthGratitude

## out of range problems -> omit offending cases
sapply(yg[,  4:9 ], function(x) sum(x < 1 | x > 7))
sapply(yg[, 10:12], function(x) sum(x < 1 | x > 5))
sapply(yg[, 13:28], function(x) sum(x < 1 | x > 9))
yg <- yg[yg$ao_4 <= 9, ]

## non-integer scale values -> round (for now)
sapply(yg[, 4:28], function(x) sum((x - round(x)) > 0))
yg[, 4:28] <- round(yg[, 4:28])


#################
## Replication ##
#################

## losd_1 (= GRAT Item 2) excluded
## gq6_6 excluded (see bottom p. 313 to top p. 314)

## Approximate replication of one of their models; could
## continue this for all of Table 1.  The test statistics are
## generally off, but the loadings are pretty close. Note that,
## when they do the measurement invariance tests following Table 1, they 
## set "std.lv=FALSE".
m1 <- cfa('f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = subset(yg, agegroup == "12-13"), std.lv = TRUE, std.ov = TRUE)
summary(m1)


##################################
## Measurement invariance tests ##
##################################

## losd_1 (= GRAT Item 2) excluded
## gq6_6 excluded (see bottom p. 313 to top p. 314)

## one-group models for GQ-6, GAC, and GRAT
m_gq6  <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, meanstructure = TRUE)
m_gac  <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3',
  data = yg, meanstructure = TRUE)
m_grat <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
 data = yg, meanstructure = TRUE)

## set up categorical and ordinal efpFunctional objects
if(!file.exists("efp-age.rda")) {
cat_age      <- catL2BB(yg$age)
cat_agegroup <- catL2BB(yg$agegroup)
set.seed(1090)
ord_age      <- ordL2BB(yg$age,      nproc = 1:4, nobs = 5000)
ord_agegroup <- ordL2BB(yg$agegroup, nproc = 1:4, nobs = 5000)
save(cat_age, cat_agegroup, ord_age, ord_agegroup, file = "efp-age.rda")
} else {
load("efp-age.rda")
}

if(FALSE) {
gefp_2 <- gefp(m2, fit = NULL, scores=estfun.lavaan,
               order.by = as.numeric(yg$agegroup), 
               vcov = info_full, sandwich = FALSE, parm=1:4)
sctest(gefp_2, functional=catL2BB(gefp_2))  # p = .01
sctest(gefp_2, functional=ordL2BB(gefp_2))  # p = .1

## Fit a one-group model for GAC and test for invariance
## in loadings:
m3 <- cfa('f1 =~ gac_1 + gac_2 + gac_3',
          data=yg, meanstructure=TRUE)

gefp_3 <- gefp(m3, fit = NULL, scores=estfun.lavaan,
               order.by = as.numeric(yg$agegroup), 
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
           f3 ~ 0*1', data=yg,
           meanstructure=TRUE)

gefp_m4 <- gefp(m4, fit = NULL, scores=estfun.lavaan,
                order.by = as.numeric(yg$agegroup), 
                vcov = info_full, sandwich = FALSE, parm=1:4)
## The catL2BB p-value is sometimes non-significant, depending on
## how one handles the errors in the data.  Also might get something
## by changing the identification constraint (choose a different
## factor loading to fix at one).
sctest(gefp_m4, functional=catL2BB(gefp_m4))  # p = .047
ordfun <- ordL2BB(gefp_m4)
sctest(gefp_m4, functional=ordfun)  # p = .001
## The plot has five categories, but there are actually six defined by
## yg$agegroup.  I haven't looked at the code that defines the ordinal
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
}
