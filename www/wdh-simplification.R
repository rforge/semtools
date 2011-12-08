## Package and data
library("lavaan")
library(strucchange)
source("lavaan-S3.R")
data("StereotypeThreat", package = "psychotools")

## Include group variable and order correspondingly
StereotypeThreat <- transform(StereotypeThreat, group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]

## Disregard groups:
m2.simp  <- c(
'ability   =~ label("load_n") * numerical  + 1 * abstract + label("load_v") * verbal',
'ability   ~  0 * 1',
'ability   ~~ label("lvar") * ability',

'numerical ~  label("mean_n:min_t") * 1',
'abstract  ~  label("mean_a:min_t") * 1',
'verbal    ~  label("mean_v:min_t") * 1',

'numerical ~~ label("var_n:min_t") * numerical',
'abstract  ~~ label("var_a:min_t") * abstract',
'verbal    ~~ label("var_v:min_t") * verbal')

m2.simp.fit <- lavaan(m2.simp, data = StereotypeThreat, meanstructure = TRUE)

gefp_m2.simp <- gefp(m2.simp.fit, fit = NULL, order.by = StereotypeThreat$gpa, sandwich = FALSE) #, parm=c(1,4,7))  To test only parameters related to numerical variable

plot(gefp_m2.simp, functional=maxBB)
plot(gefp_m2.simp, functional=meanL2BB)
