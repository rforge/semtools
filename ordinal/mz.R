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

  ## model from Merkle & Zeileis
  rval <- 'verbal =~ x1 + x2 + x3; math =~ y1 + y2 + y3'

  data <- data[, c("x1", "x2", "x3", "y1", "y2", "y3")]
  
  ## fit model
  rval <- cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE)

  ## p-value for LRT Stat
  if(!inherits(m2, "try-error") & !inherits(m3, "try-error")) {
    lrt.p <- anova(m3,m2)[[7]][2]
  } else {
    lrt.p <- NA
  }
  
  ## store (1) fitted model object, (2) lavaan coefficient
  ## labeling and ordering
  rval <- list(
    data = data,
    model = rval,
    lrt.p = lrt.p,
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
