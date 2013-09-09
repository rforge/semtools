ordfit <- function(data, silent = TRUE, suppressWarnings = TRUE, ...)
{
 ## data checks
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "y1", "y2", "y3", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("lavaan"))

  ## model from Merkle & Zeileis (for proposed tests)
  rval <- 'verbal =~ x1 + x2 + x3; math =~ y1 + y2 + y3'

  data <- data[, c("x1", "x2", "x3", "y1", "y2", "y3")]
  
  ## fit model
  rval <- try(cfa(rval, data = data, meanstructure = TRUE, std.lv = TRUE))
  
  ## store (1) fitted model object, (2) lavaan coefficient
  ## labeling and ordering
  rval <- list(
    data = data,
    model = rval,
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
