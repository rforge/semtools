
## Data-generating process
dgp <- function(nobs = 200, diff = 3, nlevels=10, gradual=FALSE, anomaly=FALSE)
{
  ## Generates data from a factor analysis model that violates
  ## measurement invariance.  Also generates an ordinal auxiliary
  ## variable related to the invariance (loosely called "age").

  ## If gradual==TRUE, mild instability starts halfway up age
  ## and gets larger.  If gradual==FALSE, the instability starts
  ## halfway up age and stays constant.
  ## If anomaly==TRUE, instability only occurs at level (1+nlevels/2).
  ## Otherwise, instability continues up the ordinal variable.
  stopifnot(require("mvtnorm"))
  
  # ses is SE increase in parameters at last ordinal category.
  # (This is constant if gradual==FALSE but not if gradual==TRUE.)
  ses <- diff
  sampsize <- nobs
  sampsize.per.nobs <- nobs %/% nlevels

  ## Deal with cases where nobs/nlevels is not an integer
  age <- rep(NA, sampsize)
  age[1:(nlevels*sampsize.per.nobs)] <- rep(1:nlevels, each=sampsize.per.nobs)
  if (any(is.na(age))) age[is.na(age)] <- sample((1:nlevels),sum(is.na(age)))
  half.level <- ifelse(nlevels %% 2 == 1, (nlevels %/% 2)+1, nlevels/2)
  
  # Define parameter vectors/matrices:
  mu <- c(29.32,24.70,14.84,10.59,19.30,18.01)

  # Loadings:
  lambda <- matrix(0,6,2)
  lambda[,1] <- c(4.92,2.96,5.96,0,0,0)
  lambda[,2] <- c(0,0,0,3.24,4.32,7.21)

  phi <- matrix(.48,2,2)
  diag(phi) <- 1

  ## Base variances
  psi <- matrix(0,6,6)
  diag(psi) <- c(26.77,13.01,30.93,3.17,8.82,22.5)

  ## Initialize data matrix:
  datmat <- matrix(NA,sampsize,6)
  colnames(datmat) <- c("x1", "x2", "x3", "y1", "y2", "y3")

  ## Generate z and u vectors for invariant individuals:
  n.invariant <- sum(age <= half.level)
  z <- t(rmvnorm(n.invariant,rep(0,2),phi))
  u <- t(rmvnorm(n.invariant,rep(0,6),psi))

  ## Based on z and u, generate data:
  tmp.ind <- which(age <= half.level)
  for (i in 1:n.invariant){
    datmat[tmp.ind[i],] <- mu + lambda%*%z[,i] + u[,i]
  }

  ## Now do the same for invariant-violating levels
  n.vary <- nlevels - half.level
  each.vary <- ses/n.vary
  for (i in (half.level+1):nlevels){
    tmp.n <- sum(age==i)
    tmp.ind <- which(age==i)
    z <- t(rmvnorm(tmp.n,rep(0,2),phi))

    mult <- ses
    if (anomaly) mult <- ifelse(i==(half.level+1), ses, 0)

    ## SES found by fitting model to population S, multiplying
    ## observed SES by sqrt(n)
    if (gradual){
      if (anomaly) stop("gradual=TRUE and anomaly=TRUE do not work together.")
      diag(psi) <- diag(psi) + (each.vary*c(62.1,26.6,82.1,9.05,19.35,51.59))/sqrt(tmp.n)
      tmp.psi <- psi
    } else {
      tmp.psi <- psi
      diag(tmp.psi) <- diag(tmp.psi) + (mult*c(62.1,26.6,82.1,9.05,19.35,51.59))/sqrt(tmp.n)
    }

    u <- t(rmvnorm(tmp.n,rep(0,6),tmp.psi))
    for (j in 1:tmp.n){
      datmat[tmp.ind[j],] <- mu + lambda%*%z[,j] + u[,j]
    }
  }
  
  cbind(as.data.frame(datmat), age)
}

## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 5000, size = 0.05, ordfun = NULL, test = NULL, verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, length(test) * nrep), ncol = length(test))
  colnames(pval) <- test
  
  for(i in 1:nrep) {
    d <- dgp(...)

    mz <- ordfit(d)

    ## Requires estfun.lavaan
    ord_gefp  <- try(gefp(mz, fit = NULL, vcov = info.mzfit, order.by = d$age, sandwich = FALSE, parm = 7:12), silent = TRUE)

    if(!inherits(ord_gefp, "try-error")) {
      pval[i, 1] <- sctest(ord_gefp,  functional = ordfun)$p.value
      pval[i, 2] <- sctest(ord_gefp,  functional = ordwmax(ord_gefp))$p.value
      pval[i, 3] <- sctest(ord_gefp,  functional = catL2BB(ord_gefp))$p.value
      #pval[i, 4] <- sctest(ord_gefp,  functional = supLM(0.1))$p.value
      pval[i, 4] <- mz$lrt.p
      pval[i, 5] <- mz$lrt.sb
      pval[i, 6] <- mz$yb.p
      pval[i, 7] <- mz$aic
    }
  }
  rval <- colMeans(pval < size, na.rm = TRUE)
  if(verbose) print(rval)
  
  return(rval)
}

## Loop over scenarios
simulation <- function(diff = seq(0, 1.5, by = 0.25),
  nobs = c(120, 480, 960), nlevels = c(4, 8, 12), gradual = FALSE,
  anomaly = FALSE, verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs)
  nprs <- nrow(prs)
  
  test <- c("ordmax","ordwmax","catdiff","lrt","lrt.sb","yb97","aic") # had suplm here
  ntest <- length(test)

  ## Only simulate critical values if we need it.
  if (file.exists("critvals.rda")){
    load("critvals.rda")
    ## Check to make sure it has what we need; if not,
    ## rerun critvals.
    if (!all(nlevels %in% as.numeric(names(cval)))) critvals(nlevels)
  } else {
    critvals(nlevels)
  }
  cval.conds <- as.numeric(names(cval))

  do.parallel <- require("parallel")
  if (do.parallel){
    pow <- mclapply(1:nprs, function(i){
      testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                nlevels = prs$nlevels[i], gradual = gradual, test = test,
                ordfun = cval[[which(cval.conds == prs$nlevels[i])]],
                anomaly = anomaly, verbose = verbose, ...)},
                    mc.cores = round(.75*detectCores()),
                    mc.preschedule = FALSE)
    pow <- t(matrix(unlist(pow), ntest, nprs))
  } else {
    pow <- matrix(rep(NA, ntest * nprs), ncol = ntest)
    for(i in 1:nprs) {
      if(verbose) print(prs[i,])

      ## Find critical values we need
      ordfun <- cval[[which(cval.conds == prs$nlevels[i])]]

      pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i], nlevels = prs$nlevels[i], ordfun = ordfun, test = test, anomaly = anomaly, verbose = verbose, ...)
    }
  }

  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(test, each = nprs))
  rval$nobs <- factor(rval$nobs)
  rval$power <- as.vector(pow)
  return(rval)
}

## Simulate critical values
critvals <- function(nlevels = c(4, 8, 12), verbose = TRUE, ...)
{
  do.parallel <- require("parallel")
  if (do.parallel){
    cval <- mclapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000),
                     mc.cores = length(nlevels))
  } else {
    cval <- lapply(nlevels, function(x) ordL2BB(rep(1, x)/x, nobs = 5000))
  }
  names(cval) <- nlevels

  save(cval, file="critvals.rda")
}
