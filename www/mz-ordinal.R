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
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances")))
  m3 <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances","residuals","means")))
  ## Get test stats at values of m3, for Satorra-Bentler correction
  ## FIXME: User-specified starting values: manipulate
  ##        parameterEstimates(m3)$est?  Must first fix a bug in lavaan.
  tmp.pars <- parameterEstimates(m3)
  ## tmp.pars$est
  ##m3@ParTable$ustart[!is.na(m3@ParTable$ustart)] <- c(coef(m3), rep(coef(m3)[7:12], 9))
  ##m2b <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
  ##              group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), start = svals, control=list(iter.max=0)))
  ## For iter.max==0, fitMeasures() will not work.  But get value
  ## of discrepancy function via m2b@Fit@fx (may need to multiply by 2)
  
  ## Get WLS fits, for Yuan-Bentler correction
  ## NB WLS is AKA ADF.  Yuan + Bentler 1997 call it Generalized
  ##    Least Squares with a weight matrix.
  ## NB It seems that these models do not converge as often
  ##    at smaller n.  The simulations will gloss over this.
  m2.wls <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances"), estimator="WLS"))

  m3.wls <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE,
            group = "age", group.equal = c("loadings","intercepts","regressions","lv.covariances","residuals","means"), estimator="WLS"))

  ## p-value for LRT Stat + Yuan-Bentler correction
  if(!inherits(m2, "try-error") & !inherits(m3, "try-error") &
     m2@Fit@converged & m3@Fit@converged) {
    lrt.p <- anova(m3,m2)[[7]][2]
  } else {
    lrt.p <- NA
  }
  if(!inherits(m2.wls, "try-error") & !inherits(m3.wls, "try-error")) {
    ## FIXME: Difference test with YB-corrected statistics may not be
    ## technically correct, but I do not see a straightforward solution.
    n <- nrow(data)

    if (m2.wls@Fit@converged){
      m2.fit <- fitMeasures(m2.wls)
      m2.yb <- m2.fit[names(m2.fit)=="chisq"]/
        (1 + m2.fit[names(m2.fit)=="chisq"]/n)
    } else {
      m2.yb <- NA
    }

    if (m3.wls@Fit@converged){
      m3.fit <- fitMeasures(m3.wls)
      m3.yb <- m3.fit[names(m3.fit)=="chisq"]/
        (1 + m3.fit[names(m3.fit)=="chisq"]/n)
    } else {
      m3.yb <- NA
    }

    yb.diff <- m3.yb - m2.yb
    yb.p <- pchisq(yb.diff, df=anova(m3,m2)[[6]][2], lower.tail=FALSE)
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
