catL2BB <- function(freq)
{
  if(inherits(freq, "gefp")) freq <- as.factor(time(freq)[-1L])
  if(is.factor(freq)) freq <- prop.table(table(freq))
  freq <- freq / sum(freq)
  ncat <- length(freq)
  
  catdiffL2 <- function(x) {
    d <- diff(c(0, x[round(cumsum(freq) * length(x))]))
    sum(d^2 / freq)
  }
  
  computePval <- function(x, nproc = 1) pchisq(x, df = (ncat - 1) * nproc, lower.tail = FALSE)
  computeCritval <- function(alpha, nproc = 1) qchisq(alpha, df = (ncat - 1) * nproc, lower.tail = FALSE)

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(time = catdiffL2, comp = sum),
    computePval = computePval,
    computeCritval = computeCritval)
}

ordL2BB <- function(freq, nobs = NULL, nproc = NULL, ...)
{
  if(inherits(freq, "gefp")) {
    if(is.null(nproc)) nproc <- NCOL(freq$process)
    freq <- prop.table(table(time(freq)[-1L]))
  } else if(is.factor(freq)) {
    if(is.null(nproc)) nproc <- 1:20
    freq <- prop.table(table(freq))
  } else {
    if(is.null(nproc)) nproc <- 1:20
  }
  freq <- freq / sum(freq)
  ncat <- length(freq)
  tcat <- cumsum(freq[-ncat])
  if(is.null(nobs)) nobs <- ceiling(50 / min(freq))

  catwmax <- function(x) {
    n <- length(x)
    tt <- seq(along = x)/n
    ix <- round(tcat * n)
    x <- x[ix]
    tt <- tt[ix]
    x <- x/(tt * (1-tt))
    return(max(x))
  }

  boundary <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    type = "b", boundary = TRUE, ...)
  {
    n <- x$nobs
    ix <- round(tcat * n)
    tt <- ix/n

    if(is.numeric(boundary)) {
      cval <- boundary
      boundary <- TRUE
    } else {
      cval <- 0
      boundary <- FALSE
    }
    bound <- cval * boundary(tt)
    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    if(aggregate) {
      proc <- zoo(rowSums(as.matrix(x$process)^2)[-1][ix], time(x)[-1][ix])
      bound <- zoo(bound, time(proc))
      tt <- zoo(tt, time(proc))      
      proc <- proc/(tt * (1-tt))
      
      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc))) ##, range(bound)))
    
      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, type = "b", ...)
      abline(0, 0)
      if(boundary) lines(bound, col = 2)	    
    } else {
      if(is.null(ylim) & NCOL(x$process) < 2) ylim <- range(c(range(x$process)))
      if(is.null(ylab) & NCOL(x$process) < 2) ylab <- "Empirical fluctuation process"

      panel <- function(x, ...)
      {
        lines(x, ...)
        abline(0, 0)
      }
      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ylim = ylim, ...)
    }
  }

  rval <- efpFunctional(lim.process = "Brownian bridge",
    functional = list(comp = function(x) sum(x^2), time = catwmax),
    boundary = boundary, plotProcess = plotProcess,
    nobs = nobs, nproc = nproc, ...)
  rval$plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    type = "b", boundary = TRUE, ...)
  {
    cval <- rval$computeCritval(alpha = alpha, nproc = NCOL(x$process))
    plotProcess(x, alpha = alpha, aggregate = aggregate,
      xlab = xlab, ylab = ylab, main = main, ylim = ylim,
      type = type, boundary = cval, ...)  
  }
  return(rval)
}

ordwmax <- function(freq, algorithm = GenzBretz(), ...)
{
  stopifnot(require("mvtnorm"))
  
  if(inherits(freq, "gefp")) {
    freq <- prop.table(table(time(freq)[-1L]))
  } else if(is.factor(freq)) {
    freq <- prop.table(table(freq))
  }
  freq <- freq / sum(freq)
  ncat <- length(freq)
  tcat <- cumsum(freq[-ncat])

  catwmax <- function(x) {
    n <- length(x)
    tt <- seq(along = x)/n
    ix <- round(tcat * n)
    x <- x[ix]
    tt <- tt[ix]
    x <- x/sqrt(tt * (1-tt))
    return(max(abs(x)))
  }

  make_sigma <- function(tt) outer(tt, tt,
    function(x, y) sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y))))))
  sigma <- make_sigma(tcat)
  
  computePval <- function(x, nproc = 1)
    (1 - pmvnorm(lower = -x, upper = x, mean = rep(0, ncat - 1), sigma = sigma)^nproc)

  computeCritval <- function(alpha, nproc = 1)
    qmvnorm((1 - alpha)^(1/nproc), tail = "both.tails", sigma = sigma)$quantile

  boundary <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL,
    boundary = TRUE, type = "b", ...)
  {
    n <- x$nobs
    ix <- round(tcat * n)
    tt <- ix/n

    cval <- if(boundary) computeCritval(alpha = alpha, nproc = NCOL(x$process)) else 0
    bound <- cval * boundary(tt)
    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    proc <- zoo(as.matrix(x$process)[-1, , drop = FALSE][ix, , drop = FALSE], time(x)[-1][ix])
    bound <- zoo(bound, time(proc))
    tt <- zoo(tt, time(proc))	   
    proc <- proc/sqrt(tt * (1-tt))
    
    if(aggregate) proc <- zoo(apply(abs(proc), 1, max), time(proc))
    
    if(is.null(ylab)) ylab <- if(aggregate) "Empirical fluctuation process" else colnames(x$process)
    if(is.null(ylim)) ylim <- range(c(range(proc), 0, cval, if(aggregate) NULL else -cval))
    
    mypanel <- function(x, ...)
    {
      lines(x, ...)
      abline(0, 0)
      if(boundary) lines(bound, col = 2)	  
      if(!aggregate & boundary) lines(-bound, col = 2)	  
    }

    plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, type = "b",
      panel = mypanel, ...)
  }

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(time = catwmax, comp = max),
    boundary = boundary,
    computePval = computePval,
    computeCritval = computeCritval,
    plotProcess = plotProcess, ...)
}

if(FALSE) {

## artificial data
set.seed(1)
d <- data.frame(z = rep(1:10, each = 30), x = runif(300), err = rnorm(300))
d$zcat <- factor(d$z, ordered = TRUE)
d$y <- with(d, ifelse(z <= 3, 0, 1 - 2 * x) + err)

## empirical fluctuation process
## (note: currently order.by still needs to be numeric, will fix in strucchange)
scus <- gefp(y ~ x, order.by = ~ z, data = d)

## supLM test (treating order.by as totally ordered)
sctest(scus, functional = supLM(0.1))

## Chi-squared test (treating ordery.by as unordered categories)
sctest(scus, functional = catL2BB(scus))

## supLM-type test (treating order.by as ordered categories)
set.seed(1)
ordfun <- ordL2BB(scus)
sctest(scus, functional = ordfun)

## corresponding visualization
plot(scus, functional = ordfun)

## comparison with supLM(0.1) test with broken ties
d$zcont <- d$z + rep(-29:0/116, 10)
scus2 <- gefp(y ~ x, order.by = ~ zcont, data = d)
plot(scus2, functional = supLM(0.1),
  col = "gray", lwd = 2, boundary = FALSE,
  xlim = c(0, 10), ylim = c(0, 20), xlab = "z")
par(new = TRUE)
plot(scus, functional = ordfun, boundary = FALSE,
  xlim = c(0, 10), ylim = c(0, 20), xlab = "", ylab = "")
abline(h = ordfun$computeCritval(alpha = 0.05, nproc = 2), col = 2)
abline(h = supLM(0.1)$computeCritval(alpha = 0.05, nproc = 2), col = 2, lty = 2)

## new (weighted) double-maximum test for ordered categories
plot(scus, functional = ordwmax(scus), aggregate = FALSE)
sctest(scus, functional = ordwmax(scus))

## example when only one parameter (intercept) changes
set.seed(1)
d$x2 <- rep(c(-1, 1), 150)
d$y2 <- with(d, ifelse(z <= 3, 0, -0.33) + err)
scus2 <- gefp(y2 ~ x2, order.by = ~ z, data = d)
sctest(scus2, functional = supLM(0.1))
sctest(scus2, functional = ordfun)
sctest(scus2, functional = ordwmax(scus2))
}
