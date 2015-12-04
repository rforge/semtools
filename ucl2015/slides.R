### R code from vignette source 'slides.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE, width = 70)
set.seed(1090)


###################################################
### code chunk number 2: packages
###################################################
library("lavaan")
library("strucchange")


###################################################
### code chunk number 3: data
###################################################
data("YouthGratitude", package = "psychotools")
compcases <- apply(YouthGratitude[, 4:28], 1,
  function(x) all(x %in% 1:9))
yg <- YouthGratitude[compcases, ]


###################################################
### code chunk number 4: models
###################################################
gq6cfa <- cfa("f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5",
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")


###################################################
### code chunk number 5: wdmo
###################################################
sctest(gq6cfa, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = "WDMo", plot = TRUE)


###################################################
### code chunk number 6: maxlmo
###################################################
sctest(gq6cfa, order.by = yg$agegroup, parm = 1:4,
  vcov = "info", functional = "maxLMo", plot = TRUE)


###################################################
### code chunk number 7: packages
###################################################
library("psychotools")
library("psychotree")


###################################################
### code chunk number 8: data
###################################################
load("MathExam.rda")
mex <- subset(MathExam, nsolved > 0 & nsolved < 13 & group == 1)


###################################################
### code chunk number 9: classes
###################################################
sapply(mex, function(x) class(x)[1])


###################################################
### code chunk number 10: solved-plot (eval = FALSE)
###################################################
## plot(mex$solved)


###################################################
### code chunk number 11: solved-plot-fig
###################################################
plot(mex$solved)


###################################################
### code chunk number 12: raschmodel-plot (eval = FALSE)
###################################################
## ram <- raschmodel(mex$solved)
## plot(ram, type = "profile")


###################################################
### code chunk number 13: raschmodel-plot-fig
###################################################
ram <- raschmodel(mex$solved)
plot(ram, type = "profile")


###################################################
### code chunk number 14: raschmodel-piplot (eval = FALSE)
###################################################
## plot(ram, type = "piplot")


###################################################
### code chunk number 15: raschmodel-piplot-fig
###################################################
plot(ram, type = "piplot")


###################################################
### code chunk number 16: raschtree (eval = FALSE)
###################################################
## rt <- raschtree(solved ~ tests + attempt + nsolved + study, data = mex)


###################################################
### code chunk number 17: raschtree-fig
###################################################
plot(rt)


