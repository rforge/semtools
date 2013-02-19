ordfit <- function(data, silent = TRUE, suppressWarnings = TRUE, ...)
{
  ## data checks
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "y1", "y2", "y3", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  ## Fit multi-group models for LRT; similar to catL2BB
  rval <- 'verbal =~ x1 + x2 + x3; math =~ y1 + y2 + y3;
           verbal ~ 0*1; math ~ 0*1'

  m2 <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), test="Satorra.Bentler"))
  m3 <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances","residuals","means"), test="Satorra.Bentler"))

  ## Get test stats at values of m3, for Satorra-Bentler 2010 correction
  ## for difference tests.  This is currently difficult in lavaan due to
  ## 1. User-defined starting values are not working correctly.  However,
  ##    you can trick lavaan by supplying a modified ParTable (filling in the
  ##    ustart column).
  ## 2. Lavaan will not return the usual S-B scaled statistic for models
  ##    that have not converged.
  ## Set up an m2 ParTable:
  ## m2b <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
  ##                group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), do.fit=FALSE))
  ## Fill this ParTable in with values from m3:
  ## m2b@ParTable$ustart <- coef(m3,"all")
  ## Supply starting values from this partable, with do.fit=FALSE
  ## m2b <- try(cfa(m2b@ParTable, data = data, meanstructure = TRUE,
  ##               std.lv = TRUE, group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), do.fit=FALSE))
  ## Problem: Still can't get a S-B scaled statistic from m2b.

  
  ## Get WLS fits, for Yuan-Bentler correction
  ## NB WLS is AKA ADF.  GLS is WLS assuming multivariate
  ##    normality (see Browne + Arminger 95, p 189).  It seems
  ##    that we need WLS here.  However, Type I errors are off
  ##    the charts.
  ## NB It seems that these models do not converge as often
  ##    at smaller n.  The simulations will gloss over this.
  zinv <- try(get.zcov(data))
  if(!inherits(zinv, "try-error")){
    m2.wls <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
                      group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), estimator="WLS", WLS.V=zinv))

    m3.wls <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
                      group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances","residuals","means"), estimator="WLS", WLS.V=zinv))
  }else{
    m2.wls <- try(y12, silent=TRUE)
    m3.wls <- try(y12, silent=TRUE)
  }

  ## p-value for LRT Stat, AIC, Sat-Bent scaled difference
  if(!inherits(m2, "try-error") & !inherits(m3, "try-error") &
     m2@Fit@converged & m3@Fit@converged) {
    lrt <- fitMeasures(m3,"chisq") - fitMeasures(m2,"chisq")
    lrt.p <- pchisq(lrt, fitMeasures(m3,"df") - fitMeasures(m2, "df"),
                    lower.tail=FALSE)
    ## This is 1 if m3 judged better, 0 otherwise.
    ## We eventually compute power as "number of aics < .05",
    ## which will give us what we want.
    aic <- fitMeasures(m3, "aic") < fitMeasures(m2, "aic")

    ## Satorra-Bentler scaled difference test automatically
    ## computed via anova()
    lrt.sb <- anova(m2,m3,SB.classic=TRUE)[[7]][2]
  } else {
    lrt.p <- NA
    aic <- NA
    lrt.sb <- NA
  }
  ## p-value for Yuan-Bentler corrected statistic
  if(!inherits(m2.wls, "try-error") & !inherits(m3.wls, "try-error")){
     # & m2.wls@Fit@converged & m3.wls@Fit@converged) {
    ## FIXME? Difference test with YB-corrected statistics problematic?
    ## Seems not for the simulation, since the data are MVN...
    n <- nrow(data)

    if (m2.wls@Fit@converged){
      m2.fit <- fitMeasures(m2.wls, "chisq")
      m2.yb <- m2.fit/(1 + m2.fit/n)
    } else {
      m2.yb <- NA
    }

    if (m3.wls@Fit@converged){
      m3.fit <- fitMeasures(m3.wls, "chisq")
      m3.yb <- m3.fit/(1 + m3.fit/n)
    } else {
      m3.yb <- NA
    }

    yb.diff <- m3.yb - m2.yb
    yb.p <- pchisq(yb.diff, df=anova(m2,m3)[[6]][2], lower.tail=FALSE)
    if(!m2.wls@Fit@converged | !m3.wls@Fit@converged) yb.p <- NA
  } else {
    yb.p <- NA
  }

  ## model from Merkle & Zeileis (for proposed tests)
  rval <- 'verbal =~ x1 + x2 + x3; math =~ y1 + y2 + y3'

  data <- data[, c("x1", "x2", "x3", "y1", "y2", "y3")]
  
  ## fit model
  rval <- cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE)
  
  ## store (1) fitted model object, (2) lavaan coefficient
  ## labeling and ordering
  rval <- list(
    data = data,
    model = rval,
    lrt.p = lrt.p,
    lrt.sb = lrt.sb,
    aic = as.numeric(aic),
    yb.p = yb.p,
    names = c("verbal=~x1", "verbal=~x2", "verbal=~x3", "math=~y1", "math=~y2", "math=~y3", "x1~~x1", "x2~~x2", "x3~~x3",
        "y1~~y1", "y2~~y2", "y3~~y3", "verbal~~math", "x1~1", "x2~1", "x3~1", "y1~1", "y2~1", "y3~1")
            )
  
  class(rval) <- "mzfit"

  return(rval)
}

coef.mzfit <- function(object, ...) {
    unclass(getMethod("coef", "lavaan")(object$model, ...))
}

logLik.mzfit <- function(object, ...) {
  rval <- as.vector(fitMeasures(object$model, "logl"))

  structure(rval, df = length(coef(object)),
    nobs = nrow(object$data), class = "logLik")
}

vcov.mzfit <- function(object, ...) {
    vc <- getMethod("vcov", "lavaan")(object$model, ...)
    class(vc) <- "matrix"
    vc
}

bread.mzfit <- function(x, ...) vcov(x) * nrow(x$data)

info.mzfit <- function(x, ...) solve(vcov(x) * nrow(x$data))

estfun.mzfit <- function(x, ...)
{
  estfun.lavaan(x$model)
}

get.zcov <- function(data){
  age <- data$age
  data <- data[,names(data) %in% c("x1", "x2", "x3", "y1", "y2", "y3")]
  mat.size <- ncol(data) + ncol(data)*(ncol(data)+1)/2

  s.zi <- NULL
  ages <- unique(age)
  for (i in 1:length(ages)){
    zmat <- t(apply(data[age==ages[i],], 1, function(x){
      zcp <- tcrossprod(x)
      c(x, as.numeric(zcp[lower.tri(zcp, diag=TRUE)]))
    }))

    n.age <- nrow(zmat)
    ## Take inverse of the covariane matrix:
    s.zinv <- solve(((n.age-1)/n.age)*cov(zmat))

    s.zi <- c(s.zi, list(s.zinv))
  }

  s.zi
}
