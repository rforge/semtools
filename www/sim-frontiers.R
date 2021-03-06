
## Data-generating process
dgp1 <- function(nobs = 240, diff = 1.5, nlevels = 3, gradual=FALSE, anomaly=FALSE, parms="loadings")
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
  theta<-c(4.92,2.96,5.96,3.24,4.32,7.21,26.77,13.01,30.93,3.17,8.82,22.5,
           -0.48,
           29.32,24.70,14.84,10.59,19.30,18.01)
  asym<-c(7.565818,4.951599,8.644422,3.093014,4.499676,7.368059,
           62.1,26.6,82.1,9.05,19.35,51.59,
           1.024326,7.139776,4.666005,8.151785,3.696972,5.242366,
            8.630417)
    # 56.672724,24.255196,74.917170,8.264380,17.664706,47.093002,
  mu<-theta[14:19]
  # Loadings:
  lambda <- matrix(0,6,2)
  lambda[,1] <- c(theta[1:3],0,0,0)
  lambda[,2] <- c(0,0,0,theta[4:6])

  phi <- matrix(theta[13],2,2)
  diag(phi) <- 1

  ## Base variances
  psi <- matrix(0,6,6)
  diag(psi) <- theta[7:12]

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

    if (parms=="loading"){parmsvector=c(1,rep(0,18))}
    if (parms=="var"){parmsvector=c(rep(0,12),1,rep(0,6))}
    if (parms=="error"){parmsvector=c(rep(0,6),1,rep(0,5),rep(0,7))}
    if (parms=="intercept"){parmsvector=c(rep(0,13),1,rep(0,5))}
   
    mult <- ses
    if (anomaly) mult <- ifelse(i==(half.level+1), ses, 0)

    ## SES found by fitting model to population S, multiplying
    ## observed SES by sqrt(n)
    if (gradual){
      if (anomaly) stop("gradual=TRUE and anomaly=TRUE do not work together.")
    theta.tmp <- theta + (each.vary*asym*parmsvector/sqrt(tmp.n))
    # Loadings:
    tmp.lambda <- matrix(0,6,2)
    tmp.lambda[,1] <- c(theta.tmp[1:3],0,0,0)
    tmp.lambda[,2] <- c(0,0,0,theta.tmp[4:6])

    tmp.phi <- matrix(theta.tmp[13],2,2)
    diag(tmp.phi) <- 1

    # variances
    tmp.psi <- matrix(0,6,6)
    diag(tmp.psi) <- theta.tmp[7:12]

    # intercepts
    tmp.mu <- theta.tmp[14:19]
    } else {
    theta.tmp <- theta + (mult*asym*parmsvector/sqrt(tmp.n))
    # Loadings:
    tmp.lambda <- matrix(0,6,2)
    tmp.lambda[,1] <- c(theta.tmp[1:3],0,0,0)
    tmp.lambda[,2] <- c(0,0,0,theta.tmp[4:6])

    tmp.phi <- matrix(theta.tmp[13],2,2)
    diag(tmp.phi) <- 1

    # variances
    tmp.psi <- matrix(0,6,6)
    diag(tmp.psi) <- theta.tmp[7:12]  
    }

    # intercepts
    tmp.mu <- theta.tmp[14:19]
    
    u.tmp <- t(rmvnorm(tmp.n,rep(0,6),tmp.psi))
    z.tmp <- t(rmvnorm(tmp.n,rep(0,2),tmp.phi))
    for (j in 1:tmp.n){
      datmat[tmp.ind[j],] <- tmp.mu + tmp.lambda%*%z.tmp[,j] + u.tmp[,j]
    }
  }
  
  cbind(as.data.frame(datmat), age)
}

dgp2 <-function(nobs = 240, diff = 1.5, nlevels = 3, gradual=FALSE,
                anomaly=FALSE, parms="loadings")
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
  theta<-c(4.92,2.96,5.96,1,3.24,4.32,7.21,26.77,13.01,30.93,3.17,8.82,22.5,-0.48,
           29.32,24.70,14.84,10.59,19.30,18.01)
  asym<-c(10.037376,4.978146,8.868835,9.009867,3.095191,4.501605,
           7.371680,60.969079,24.444897,78.555777,8.298542,17.701359,
          47.210688,1.127309,6.874096,4.666005,8.151785,3.696972,
           5.242366,8.630417)
  mu<-theta[15:20]
  # Loadings:
  lambda <- matrix(0,6,2)
  lambda[,1] <- c(theta[1:3],0,0,0)
  lambda[,2] <- c(theta[4],0,0,theta[5:7])

  phi <- matrix(theta[14],2,2)
  diag(phi) <- 1

  ## Base variances
  psi <- matrix(0,6,6)
  diag(psi) <- theta[8:13]

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

  if(parms=="extra"){parmsvector=c(rep(0,3),1,rep(0,3),rep(0,6),0,rep(0,6))}
  if(parms=="extra+loading"){parmsvector=c(1,rep(0,2),0,rep(0,3),rep(0,6),0,rep(0,6))}
  if(parms=="extra+var"){parmsvector=c(rep(0,3),0,rep(0,3),rep(0,6),1,rep(0,6))}
  if(parms=="extra+error"){parmsvector=c(rep(0,3),0,rep(0,3),1,rep(0,5),0,
                               rep(0,6))}
  if(parms=="extra+intercept"){parmsvector=c(rep(0,14),
                                   1, rep(0,5))}
 
   
    mult <- ses
    if (anomaly) mult <- ifelse(i==(half.level+1), ses, 0)

    ## SES found by fitting model to population S, multiplying
    ## observed SES by sqrt(n)
    if (gradual){
      if (anomaly) stop("gradual=TRUE and anomaly=TRUE do not work together.")
    theta.tmp <- theta + (each.vary*asym*parmsvector/sqrt(tmp.n))
    # Loadings:
    tmp.lambda <- matrix(0,6,2)
    tmp.lambda[,1] <- c(theta.tmp[1:3],0,0,0)
    tmp.lambda[,2] <- c(theta.tmp[4],0,0,theta.tmp[5:7])

    tmp.phi <- matrix(theta.tmp[14],2,2)
    diag(tmp.phi) <- 1

    # variances
    tmp.psi <- matrix(0,6,6)
    diag(tmp.psi) <- theta.tmp[8:13]

    # intercepts:
    tmp.mu <- theta.tmp[15:20]
    } else {
    theta.tmp <- theta + (mult*asym*parmsvector/sqrt(tmp.n))
    # Loadings:
    tmp.lambda <- matrix(0,6,2)
    tmp.lambda[,1] <- c(theta.tmp[1:3],0,0,0)
    tmp.lambda[,2] <- c(theta.tmp[4],0,0,theta.tmp[5:7])

    tmp.phi <- matrix(theta.tmp[14],2,2)
    diag(tmp.phi) <- 1

    # variances
    tmp.psi <- matrix(0,6,6)
    diag(tmp.psi) <- theta.tmp[8:13]  
    }

    # intercepts:
    tmp.mu <- theta.tmp[15:20]
    
    u.tmp <- t(rmvnorm(tmp.n,rep(0,6),tmp.psi))
    z.tmp <- t(rmvnorm(tmp.n,rep(0,2),tmp.phi))
    for (j in 1:tmp.n){
      datmat[tmp.ind[j],] <- tmp.mu + tmp.lambda%*%z.tmp[,j] + u.tmp[,j]
    }
  }
  
  cbind(as.data.frame(datmat), age)
}



## Evaluate power simulation on a single dgp() scenario
testpower <- function(nrep = 5000, size = 0.05, ordfun = NULL, parnum = parnum,
                      test = test, sim="sim1", verbose = TRUE, ...)
{
  pval <- matrix(rep(NA, length(test) * nrep), ncol = length(test))
  colnames(pval) <- test
  dfun <- paste("dgp", substr(sim,4,4), sep="")
  
  for(i in 1:nrep) {
    d <- do.call(dfun, list(...))
    mz <- ordfit(d)

    ## Requires estfun.lavaan
    ord_gefp <-vector("list",length(parnum))
    for (j in 1: length(parnum)){
    ord_gefp[[j]]  <- try(gefp(mz, fit = NULL, vcov = info.mzfit,
                               order.by = d$age, sandwich =FALSE,
                               parm = eval(parse(text=parnum[j]))),
                          silent = TRUE)

    if(!inherits(ord_gefp[[j]], "try-error")) {
      pval[i, (j-1)*3+1] <- sctest(ord_gefp[[j]], functional = ordfun)$p.value
      pval[i, (j-1)*3+2] <- sctest(ord_gefp[[j]], functional = ordwmax(ord_gefp[[j]]))$p.value
      pval[i, (j-1)*3+3] <- sctest(ord_gefp[[j]], functional = catL2BB(ord_gefp[[j]]))$p.value
    }
  }
  }
  rval <- colMeans(pval < size, na.rm = TRUE)
  if(verbose) print(rval)
  
  return(rval)
}



## Loop over scenarios
simulation <- function(diff = seq(0, 4, by = 0.25),
                       parms=c("error","loading","var"),
                       nobs = c(120, 480, 960), nlevels = c(4, 8, 12),
                       gradual = FALSE, anomaly = FALSE, verbose = TRUE, ...)
{
  prs <- expand.grid(diff = diff, nlevels = nlevels, nobs = nobs, parms = parms)
  nprs <- nrow(prs)

  parnum <- c(1,13,7,14,"1:6","7:12", "14:19")

  tname <- c("ordmax","ordwmax","catdiff") # had suplm here
  
  test <- paste(rep(tname,length(parnum)), rep(parnum, each=length(tname)),sep="")
  ntest <- length(test)

  ## Only simulate critical values if we need it.
  if (file.exists("critvals.rda")){
    load("critvals.rda")
    ## Check to make sure it has what we need; if not,
    ## rerun critvals.
    if (!all(nlevels %in% as.numeric(names(cval)))) critvals(nlevels)
  } else {
    cval <- critvals(nlevels)
    save(cval, file="critvals.rda")
  }
  cval.conds <- as.numeric(names(cval))

  pow <- matrix(rep(NA,ntest*nprs),ncol=ntest)
  do.parallel <- require("parallel")

  if (do.parallel){
    pow <- mclapply(1:nprs, function(i){
      testpower(diff = prs$diff[i], nobs = prs$nobs[i], parms = prs$parms[i],
                parnum = parnum, nlevels = prs$nlevels[i], gradual = gradual,
                test = test,
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

      pow[i,] <- testpower(diff = prs$diff[i], nobs = prs$nobs[i],
                           parms = prs$parms[i],parnum = parnum,
                           nlevels = prs$nlevels[i], ordfun = ordfun,
                           test = test, anomaly = anomaly,
                           verbose = verbose, ...)
    }
  }

  rval <- data.frame()
  for(i in 1:ntest) rval <- rbind(rval, prs)
  rval$test <- factor(rep(rep(tname,length(parnum)), each = nprs),
                      levels=c("ordmax","ordwmax","catdiff"))
  rval$pars <- factor(rep(rep(parnum,each=length(tname)),each=nprs),
                      levels=parnum)
  rval$nobs <- factor(rval$nobs)
  rval$parms <- factor(rval$parms)
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



## Examples of running the simulation
if(FALSE){

library("lavaan")
library("strucchange")
library("mvtnorm")
source("sim-frontiers.R")
source("mz-frontiers.R")


## seed for replication
  RNGkind(kind = "default", normal.kind = "default")
  set.seed(1090)
  

  
## Short versions:
## Simulation 1:
simtry <- simulation(sim = "sim1", nobs = 480, nrep = 300,
                     diff = c(1, 4), nlevels = 8,
                     parms = "intercept")
## Simulation 2:
simtry <- simulation(sim = "sim2", nobs = 480, nrep=300,
                     diff = c(1, 4), nlevels = 8,
                     parms = "extra", "extra+intercept")

## Save some results to file
save(simtry,file ="simtry.rda")
 
}
