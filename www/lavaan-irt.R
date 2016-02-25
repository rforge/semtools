mzirtfit.pml <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "x4", "x5", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  data <- data[, c("x1", "x2", "x3", "x4", "x5")]

  ## fit model
  f05 <- paste(paste("lam", 1:5, "*", names(data)[1:5],
  sep = ""), collapse = " + ")
  f05 <- paste("theta =~ ", f05, sep = "")
  f06 <- paste(names(data)[1:5], " | th", 1:5, "*t1",
  sep = "")
  f06 <- paste(f06, collapse = " \n ")
  f07 <- paste("alpha", 1:5, " := (lam",1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f07 <- paste(f07, collapse = " \n ")
  f08 <- paste("beta", 1:5, " := (-th", 1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f08 <- paste(f08, collapse = " \n ")

mod <- paste(f05, f06, f07, f08, sep = " \n ")
writeLines(mod, "mod.txt")

  rval <- try(cfa(mod, data=data, std.lv=TRUE,
          ordered=names(data)[1:5], estimator="PML"))
  return(rval)
}

mzgroupirtfit.pml <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "x4", "x5", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  data <- data[, c("x1", "x2", "x3", "x4", "x5", "age")]

  ## fit model
  f05 <- paste(paste("lam", 1:5, "*", names(data)[1:5],
  sep = ""), collapse = " + ")
  f05 <- paste("theta =~ ", f05, sep = "")
  f06 <- paste(names(data)[1:5], " | th", 1:5, "*t1",
  sep = "")
  f06 <- paste(f06, collapse = " \n ")
  f07 <- paste("alpha", 1:5, " := (lam",1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f07 <- paste(f07, collapse = " \n ")
  f08 <- paste("beta", 1:5, " := (-th", 1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f08 <- paste(f08, collapse = " \n ")
  f09 <- paste(names(data)[1:5], " ~*~ 1*", names(data)[1:5],
               sep = "", collapse = " \n ")
  f10 <- "theta ~~ c(1, NA, NA, NA, NA, NA, NA, NA)*theta"
  f11 <- "theta ~ c(0, NA, NA, NA, NA, NA, NA, NA)*1"
  simmodgroup <- paste(f05, f06, f07, f08, f09, f10, f11, sep = " \n ")
  writeLines(simmodgroup, "simmodgroup.txt")

  rval <- try(cfa(simmodgroup, data=data, std.lv=TRUE,
          ordered=names(data)[1:5], group = "age", group.equal=
          c("loadings", "thresholds"), estimator="PML"))
  return(rval)
}

mzgroupmeanirtfit.pml <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "x4", "x5", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  data <- data[, c("x1", "x2", "x3", "x4", "x5", "age")]

  ## fit model
  f05 <- paste(paste("lam", 1:5, "*", names(data)[1:5],
  sep = ""), collapse = " + ")
  f05 <- paste("theta =~ ", f05, sep = "")
  f06 <- paste(names(data)[1:5], " | th", 1:5, "*t1",
  sep = "")
  f06 <- paste(f06, collapse = " \n ")
  f07 <- paste("alpha", 1:5, " := (lam",1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f07 <- paste(f07, collapse = " \n ")
  f08 <- paste("beta", 1:5, " := (-th", 1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f08 <- paste(f08, collapse = " \n ")
  f09 <- paste(names(data)[1:5], " ~*~ 1*", names(data)[1:5],
               sep = "", collapse = " \n ")
  f10 <- "theta ~ c(0, NA, NA, NA, NA, NA, NA, NA)*1"
  f11 <- "theta ~~ c(1, 1, 1, 1, 1, 1, 1, 1)*theta"
  simmodgroupmean <- paste(f05, f06, f07, f08, f09, f10, f11, sep = " \n ")
  writeLines(simmodgroupmean, "simmodgroupmean.txt")

  rval <- try(cfa(simmodgroupmean, data=data, std.lv=TRUE,
          ordered=names(data)[1:5], group = "age", group.equal=
          c("loadings", "thresholds"), estimator="PML"))
  return(rval)
}

mzgrconsirtfit.pml <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "x4", "x5", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  data <- data[, c("x1", "x2", "x3", "x4", "x5", "age")]

  ## fit model
  f05 <- paste(paste("lam", 1:5, "*", names(data)[1:5],
  sep = ""), collapse = " + ")
  f05 <- paste("theta =~ ", f05, sep = "")
  f06 <- paste(names(data)[1:5], " | th", 1:5, "*t1",
  sep = "")
  f06 <- paste(f06, collapse = " \n ")
  f07 <- paste("alpha", 1:5, " := (lam",1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f07 <- paste(f07, collapse = " \n ")
  f08 <- paste("beta", 1:5, " := (-th", 1:5, ")/sqrt(1-lam",
  1:5, "^2)", sep = "")
  f08 <- paste(f08, collapse = " \n ")
  f09 <- paste(names(data)[1:5], " ~*~ 1*", names(data)[1:5],
               sep = "", collapse = " \n ")
  f10 <- "theta ~~ c(1, NA, NA, NA, NA, NA, NA, NA)*theta +
                   c(NA, 'v', 'v', 'v', 'v', 'v', 'v', 'v')*theta "
  f11 <- "theta ~ c(0, NA, NA, NA, NA, NA, NA, NA)*1 +
                  c(NA, 'm', 'm', 'm', 'm', 'm', 'm', 'm')*1"
  simmodgrcons <- paste(f05, f06, f07, f08, f09, f11, f10, sep = " \n ")
  writeLines(simmodgrcons, "simmodgrcons.txt")

  rval <- try(cfa(simmodgrcons, data=data, std.lv=TRUE,
          ordered=names(data)[1:5], group = "age", group.equal=
          c("loadings", "thresholds"), estimator="PML"))
  return(rval)
}

estfun.lavaan <- pmlScores <- function(object){
  library("lavaan")
  ## pull arguments value from objects
  ## shortcuts
  lavdata <- object@Data
  lavmodel <- object@Model
  lavsamplestats <- object@SampleStats
  lavfit <- object@Fit
  lavcache <- object@Cache

  ## number variables/sample size
  ntab <- unlist(lavsamplestats@nobs)
  ntot <- lavsamplestats@ntotal

  Score.mat <- matrix(NA, ntot, length(coef(object)))

  for (g in 1:object@SampleStats@ngroups){
    score.H1 <- lavaan:::pml_deriv1(Sigma.hat=
                         lavfit@Sigma.hat[[g]],
                         TH=lavfit@TH[[g]],
                         th.idx=lavmodel@th.idx[[g]],
                         num.idx=lavmodel@num.idx[[g]],
                         X=lavdata@X[[g]],
                         eXo=lavdata@eXo[[g]],
                         lavcache = lavcache[[g]],
                         scores=TRUE, negative = TRUE)
  ## Transform the score matrix w.r.t. estimated parameters
  ## instead of correlations.
  ## Delta is the matrix dealing with chain rule.
    Delta <- lavaan:::computeDelta(lavmodel)[[g]]
    wi <- lavdata@case.idx[[g]]
    Score.mat[wi,] <- -score.H1 %*% Delta
    colnames(Score.mat) <- names(coef(object))
  }
   if(nrow(object@Model@eq.constraints.K)==0){
   Score.mat} else{
   Score.mat <- Score.mat
   #Score.mat %*% object@Model@eq.constraints.K
   }
   Score.mat
}


## Hessian is not equal to the crossproduct of scores. We need to get the crossproduct of scores
## directly from frist order function B0.group[[g]].

info_full <- function(object, ...){
    vcov <-lavaan:::lav_model_vcov #robust.huber.white
    varcov <- vcov(lavmodel = object@Model,
                   lavsamplestats = object@SampleStats,
                   lavoptions = object@Options,
                   lavdata = object@Data,
                   lavpartable = object@ParTable,
                   lavcache = object@Cache,control=list())
    attr(varcov, "B0.group")[[1]]
}

