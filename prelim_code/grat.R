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

## non-integer scale values -> omit offending cases
## NB: Critical values need to be re-simulated
## (Could also make them NA, then use missing="listwise"
##  cfa() commands to exclude cases with missing values.
##  However, this gives us a different sample size for
##  each model, meaning we have to re-simulate the critical
##  value for each model.)
##non.int <- (yg[,4:28] - round(yg[,4:28])) > 0
##n.non.int <- apply(non.int, 1, sum)
##yg <- yg[n.non.int==0,]


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

## one-group models for GAC, GRAT, and GQ-6
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
m_gq6  <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, meanstructure = TRUE)

## set up categorical and ordinal efpFunctional objects
## for age groups 10-11, 12, 13, 14, 15, 16, 17-19
if(!file.exists("efp-age.rda")) {
  ## relative frequencies
  freq_age <- prop.table(table(
    cut(yg$age, breaks = c(-Inf, 11:16 + 0.5, Inf), labels = 11:17)))

  ## supLM test
  sup_age <- supLM(freq_age[1], 1 - freq_age[7])

  ## Chi-squared test
  cat_age <- catL2BB(freq_age)

  ## supLM test for ordered categories
  set.seed(1090)
  ord_age <- ordL2BB(freq_age, nproc = 1:6, nobs = 5000)

  save(freq_age, sup_age, cat_age, ord_age, file = "efp-age.rda")
} else {
  load("efp-age.rda")
}

## measurement invariance fluctuation tests
mitests <- function(object, parm = NULL, plot = TRUE) {
  ## fluctuation process
  gefp_age <- gefp(object, fit = NULL, order.by = yg$age, 
    vcov = info_full, sandwich = FALSE, parm = parm)

  ## p-values
  rval <- c(
    sctest(gefp_age,	  functional = supLM(0.1))$p.value,
    sctest(gefp_age,	  functional = cat_age)$p.value,
    sctest(gefp_age,	  functional = ord_age)$p.value)
  names(rval) <- c("supLM", "catL2BB", "ordL2bb")
  
  ## visualizations
  if(plot) {
    par(mfrow = c(1, 2))
    plot(gefp_age, functional = cat_age, xlab = "Parameter", main = "Chi-squared test")
    stat  <- ord_age$computeStatistic(gefp_age$process)
    cval  <- ord_age$computeCritval(alpha = 0.05, nproc = length(parm))
    cval2 <- sup_age$computeCritval(alpha = 0.05, nproc = length(parm))
    plot(gefp_age, functional = ord_age, ylim = c(0, max(cval, cval2, stat)),
      xlab = "Age", main = "supLM test (discrete/continuous)")
    abline(h = cval,  col = 2)
    abline(h = cval2, col = 2, lty = 2)
  }
  
  return(rval)
}

## GAC factor loadings
mitests(m_gac,  parm = 1:2)
## -> all non-significant
## -> p-values in the expected order:
## unordered categorical > ordered numeric > ordered categorial

## GRAT factor loadings
mitests(m_grat, parm = 1:4)
mitests(m_grat, parm = 5:9)
mitests(m_grat, parm = 10:12)
## -> LOSD subscale seems to violate measurement invariance
## -> SA and AO appear to be ok

## GQ-6 factor loadings
mitests(m_gq6,  parm = 1:4)
## -> somewhat surprising results because ordered categorical
## is non-significant while the other two are significant


####################
## GQ-6 in detail ##
####################

## Question 1: Why is supLM significant but ordL2BB not?
## Answer: Coincidence due to ordering within the age groups.
## Illustration: Break ties randomly and compute average p-value.
set.seed(1)
p_gq6 <- sapply(1:100, function(iii) sctest(gefp(m_gq6, fit = NULL,
  order.by = yg$age + runif(nrow(yg), -0.1, 0),
  vcov = info_full, sandwich = FALSE, parm = 1:4),
  fun = sup_age)$p.value)
summary(p_gq6)
## -> the p-values vary a lot but are non-significant "on average"
## if the exact birth-date were known one could have a closer look...


## Question 2: Why is the unordered categorical test clearly significant?
## Is this due to a specific item loading or specific age groups?
## Answer: The effect does not seem to be monotonic in age. The model has
## a relatively poor for f1=~gq6_3 at ages 12, 15, 17 and for
## f1=~gq6_2 and f1=~gq6_4 at age 16. Hence, the ordered tests have trouble
## picking this up.
## Illustration: Look at the contributions to the chi-squared statistic.

## empirical fluctuation process and categorical chi-squared statistic
scus_gq6 <- gefp(m_gq6, fit = NULL, order.by = yg$age,
  vcov = info_full, sandwich = FALSE, parm = 1:4)
sctest(scus_gq6, functional = cat_age)

## compute contributions to statistic "by hand"
## disaggregate cumulative sums
s_gq6 <- diff(coredata(scus_gq6$process))
## compute sums within age groups (merging <= 11 and >= 17)
s_gq6 <- aggregate(s_gq6,
  list(pmax(11, pmin(17, time(scus_gq6$process)[-1]))), sum)
rownames(s_gq6) <- s_gq6[, 1]
s_gq6 <- s_gq6[, -1]
## the contributions are weighted squares of the sums
s_gq6 <- s_gq6^2 / as.vector(freq_age)
## replicate test statistic and p-value of Chi-squared test
sum(s_gq6)
pchisq(sum(s_gq6), (7 - 1) * 4, lower.tail = FALSE)

## the contributions > 3.84 = qchisq(0.95, 1) are the five cells mentioned above
round(s_gq6, digits = 1)

#############################################
## Try a new set of parameter restrictions ##
#############################################

## one-group models for GAC, GRAT, and GQ-6
## (same as before, but order of osberved variables
##  is changed in the equations)
m_gac2  <- cfa(
  'f1 =~ gac_2 + gac_1 + gac_3',
  data = yg, meanstructure = TRUE)
m_grat2 <- cfa(
  'f1 =~ losd_3 + losd_2 + losd_4 + losd_5 + losd_6
   f2 =~ sa_2 + sa_1 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_2 + ao_1 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE)
m_gq62  <- cfa(
  'f1 =~ gq6_2 + gq6_1 + gq6_3 + gq6_4 + gq6_5',
  data = yg, meanstructure = TRUE)

## GAC factor loadings
mitests(m_gac2,  parm = 1:2)
## -> all non-significant
## -> p-values are now:
## unordered categorical > ordered categorical > ordered numeric

## GRAT factor loadings
mitests(m_grat2, parm = 1:4)
mitests(m_grat2, parm = 5:9)
mitests(m_grat2, parm = 10:12)
## -> LOSD subscale more clearly violates measurement invariance
## -> SA and AO still ok

## GQ-6 factor loadings
mitests(m_gq62,  parm = 1:4)
## -> ordered categorical is now just significant, and
## the other two are also significant.


###########################################
## Try other restrictions, focus on GRAT ##
###########################################
m_grat3 <- cfa(
  'f1 =~ losd_4 + losd_2 + losd_3 + losd_5 + losd_6
   f2 =~ sa_2 + sa_1 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_2 + ao_1 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE)

m_grat4 <- cfa(
  'f1 =~ losd_5 + losd_2 + losd_4 + losd_3 + losd_6
   f2 =~ sa_2 + sa_1 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_2 + ao_1 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE)

m_grat5 <- cfa(
  'f1 =~ losd_6 + losd_2 + losd_4 + losd_5 + losd_3
   f2 =~ sa_2 + sa_1 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_2 + ao_1 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE)

## Restrict latent variance to be 1
m_grat6 <- cfa(
  'f1 =~ losd_6 + losd_2 + losd_4 + losd_5 + losd_3
   f2 =~ sa_2 + sa_1 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_2 + ao_1 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, std.lv=TRUE, meanstructure = TRUE)

mitests(m_grat3, parm = 1:4)
mitests(m_grat4, parm = 1:4)
mitests(m_grat5, parm = 1:4)
## An extra loading is estimated with this parameterization:
mitests(m_grat6, parm = 1:5)

###########################################
## Testing within multiple group models ###
###########################################
## Define the grouping variable
yg <- transform(yg, age2=cut(age, breaks = c(-Inf, 11:16 + 0.5, Inf), labels = 11:17))

m_grat_grp <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE, group="age2",
  group.equal="loadings")

## 1:4 for losd; 5:9 for sa; 10:12 for ao
## losd loadings: ordinal test has smaller p-value,
##      but none approach significance
mitests(m_grat_grp, parm=1:4)
## sa loadings: nothing
mitests(m_grat_grp, parm=5:9)
## ao loadings: nothing
mitests(m_grat_grp, parm=10:12)

## Now refit model with residuals equal as well
m_grat_grp2 <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4
   f1 ~ 0*1
   f2 ~ 0*1
   f3 ~ 0*1',
  data = yg, meanstructure = TRUE, group="age2",
  group.equal=c("loadings","residuals"))

## losd residuals: all significant
mitests(m_grat_grp2, parm=13:17)
## sa residuals: critical value hasn't been
##    simulated (the code above has already been modified,
##    need to re-simulate).
mitests(m_grat_grp2, parm=18:23)
## ao residuals: ordinal and supLM not significant, cat is
mitests(m_grat_grp2, parm=24:27)

####################################################
## Now try a similar procedure on the other measures

## gac:
m_gac_grp  <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3',
  data = yg, meanstructure = TRUE, group="age2",
  group.equal="loadings")

## Loadings are significant with ordinal test,
## not with others
mitests(m_gac_grp, parm=1:2)


## gq6:
m_gq6_grp  <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, meanstructure = TRUE, group="age2",
  group.equal="loadings")

## All are significant:
mitests(m_gq6_grp, parm=1:4)
