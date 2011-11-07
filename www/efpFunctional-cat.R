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
  if(is.null(nobs)) nobs <- 50 / min(freq)

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
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL, type = "b", ...)
  {
    n <- x$nobs
    ix <- round(tcat * n)
    tt <- ix/n

##    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary(tt)
    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    if(aggregate) {
      proc <- zoo(rowSums(as.matrix(x$process)^2)[-1][ix], time(x)[-1][ix])
##      bound <- zoo(bound, time(proc))
      tt <- zoo(tt, time(proc))      
      proc <- proc/(tt * (1-tt))
      
      if(is.null(ylab)) ylab <- "Empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc))) ##, range(bound)))
    
      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, type = "b", ...)
      abline(0, 0)
##      lines(bound, col = 2)	    
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

  efpFunctional(lim.process = "Brownian bridge",
    functional = list(comp = function(x) sum(x^2), time = catwmax),
    boundary = boundary, plotProcess = plotProcess,
    nobs = nobs, nproc = nproc, ...)
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
abline(h = ordfun$computeCritval(alpha = 0.05, nproc = 2), col = 2)

}
