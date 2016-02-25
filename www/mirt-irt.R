## mml using ltm package
mzirtfit.ltm <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
  data <- as.data.frame(data)
  stopifnot(c("x1", "x2", "x3", "x4", "x5", "age") %in% names(data))
  age <- data$age

  ## require package
  stopifnot(require("ltm"))

  data <- data[, c("x1", "x2", "x3", "x4", "x5")]
  
  ## fit model
  rval <- try(ltm(data ~ z1))
  return(rval)
}


estfun.ltm <- ltmScores <- function(object)
{ 
  tm <- terms(object$formula)
  av <- all.vars(object$formula)
  factors <- sum(av %in% c("z1", "z2"))
  if (factors > 2)
        stop("\nMaximum number of factors to include is 2.")
  if ((factors == 1 & !"z1" %in% av) || (factors == 2 & !c("z1", "z2") %in% av))  
        stop("\nyou have to use 'z1' for the first factor and 'z2' for the second one.")
  tm.lab <- attr(tm, "term.labels")
  ltst <- list(factors = factors, inter = "z1:z2" %in% tm.lab, quad.z1 = "I(z1^2)" %in% tm.lab, 
                 quad.z2 = "I(z2^2)" %in% tm.lab)
  betas <- object$coefficient
  p <- nrow(object$coefficient)
  q. <- 1 + factors + sum(unlist(ltst[-1]))
  X <- as.matrix(object$X)
  mX <- 1-X
  GHw <- object$GH[[2]]
  Z <- object$GH[[1]]
  obs <- nrow(X)
  if (any(na.ind <- is.na(X)))
  X[na.ind] <- mX[na.ind] <- 0
  constraint <- object$constraint

  library("ltm")
  betas.ltm <- ltm:::betas.ltm
  probs <- ltm:::probs

  pr <- probs(Z %*% t(betas))
  p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
  p.x <- c(p.xz %*% GHw)
  p.zx <- (p.xz / p.x) * obs
  #scores <- matrix(0, p, q.)
  scores <- vector("list", p)
    for (i in 1:p) {
      ind. <- na.ind[, i]
      Y <- outer(X[, i], pr[, i], "-")
      Y[ind., ] <- 0
      #scores[i, ] <- -colSums((p.zx * Y) %*% (Z * GHw))
      scores[[i]] <- -(p.zx * Y) %*% (Z * GHw)
    }
      #create matrix contains scores for each individual(n), each parameter
      #(p*q.)
      scoresuse <- matrix(unlist(scores), obs, p*q.)
    if (!is.null(constraint))
        scoresuse[-((constraint[, 2] - 1) * p + constraint[, 1])]
        # return score matrix for each individual, each parameter.
        # permutate the scoresuse matrix to match up with the parms order in
        # vcov function.
    ## will mess up with the constraint. Don't do constraint.
    scorespermu <- matrix(NA, obs,p*q.)
    for (i in 1:p){
           scorespermu[,i] <- scoresuse[,(2*i-1)]
           scorespermu[,p+i] <- scoresuse[,(2*i)]
       }
  scorespermu
}

loglikltm_i <-
function (object, betas, i) {
  # conditions
  # feed in the value needed.
  p <- nrow(object$coefficient)
  q.<- 2
  Z <- object$GH[[1]]
  X <- as.matrix(object$X)
  mX <- 1-X
  GHw <- object$GH[[2]]
  obs <- nrow(X)
  constraint <- NULL

  # functions
    probs <- ltm:::probs
    pr <- probs(Z %*% t(betas))
    # generate p.xz for each row instead of every observation in X. This is
    # the only change from original loglikltm function.
    p.xz <- exp(X[i,] %*% t(log(pr)) + mX[i,] %*% t(log(1-pr)))
    p.x <- rep(c(p.xz %*% GHw),obs)
    scorefun <- -sum(log(p.x))
    return(scorefun)
}

# in lavaan vcov is (1/N)*cov, therefor info_full=solve(1/N*cov*N)
# in ltm vcov is hessian, which is N*cov, therefore info_full=solve(N*cov/N)
info_full <- function(object,...)solve(vcov(object)/nrow(object$X))

## mml using mirt package. 
mzirtfit.mirt <- function(data, silent = TRUE, suppressWarnings = TRUE, ...){
    data$grcons <- as.factor(ifelse(as.numeric(data$age) > 1,1,0))
    models <- 'F1 = 1-5'
    mmlrval<- try(multipleGroup(data[,1:5], models, SE=FALSE,
                                group = data$grcons,
                                invariance=c('slopes',
                                'intercepts','free_var','free_means'), 
                                itemtype='2PL'))
    return(mmlrval)
}

## define functions
jointdist <- function(theta, x, alph, bet, mu, std){
    res <- matrix(NA, length(theta), length(alph))
    for(j in 1:length(theta)){
        res[j,] <- dbinom(x, size=1, prob=plogis(alph*theta[j] + bet))
    }
    apply(res, 1, prod) * dnorm(theta, mu, std)
}

# deriv wrt slope
fun1 <- function(theta, x, alph, bet, xvec, avec, bvec, mu, std){
    (theta * (x - plogis(alph*theta + bet))) *
        jointdist(theta, xvec, avec, bvec, mu, std)
}

## deriv wrt intercept
fun2 <- function(theta, x, alph, bet, xvec, avec, bvec, mu, std){
    (x - plogis(alph*theta + bet)) * jointdist(theta, xvec, avec, bvec, mu, std)
}

## deriv wrt mu
fun3 <- function(theta, x, alph, bet, xvec, avec, bvec, mu, std){
    ((theta-mu)*std^(-2)) * jointdist(theta, xvec, avec, bvec, mu, std)
}

## deriv wrt sigma
fun4 <- function(theta, x, alph, bet, xvec, avec, bvec, mu, std){
    (0.5*(theta-mu)^2*std^(-4)-0.5*std^(-2)) * jointdist(theta, xvec, avec, bvec, mu, std)
}

## get score matrix
estfun.mirt <- function(mmlrval){
  data <- mmlrval@Data$data
  totalpar <- length(as.numeric(sapply(coef(mmlrval)[[2]], function(x) x[1:2])))
  itempars <- as.numeric(sapply(coef(mmlrval)[[2]], function(x) x[1:2]))[1:(2*ncol(data))]
  personpars <- as.numeric(sapply(coef(mmlrval)[[2]], function(x) x[1:2]))[(2*ncol(data)+1):totalpar]
 
  score2 <- matrix(NA, nrow(data), totalpar)
  ## ML estimates of item parameters:

  avec <- itempars[((1:ncol(data))*2 - 1)]  # slopes
  bvec <- itempars[(1:ncol(data))*2] # intercepts

  ## baseline group
  for(i in 1:nrow(data[which(mmlrval@Data$group==0),])){
    normconst <- integrate(jointdist, -6, 6, x= as.numeric(data[i,1:ncol(data)]), 
                           alph=avec, bet=bvec, mu=0, std=1)$value

    for(j in 1:ncol(data)){
      par1loc <- (j-1)*2 + 1

      score2[i,par1loc] <- (1/normconst) * integrate(fun1, -6, 6, x=data[i,j],
                                             alph=itempars[par1loc],
                                             bet=itempars[(par1loc+1)],
                                             xvec = as.numeric(data[i,1:ncol(data)]), 
                                             avec=avec,
                                             bvec=bvec, mu=0, std=1)$value

      score2[i,(par1loc+1)] <- (1/normconst) * 
                                   integrate(fun2, -6, 6, x=data[i,j],
                                             alph=itempars[par1loc],
                                             bet=itempars[(par1loc+1)],
                                             xvec = as.numeric(data[i,1:ncol(data)]), 
                                             avec=avec, 
                                             bvec=bvec, mu=0, std=1)$value
    }
    score2[i,(totalpar-1)] <- 0
    score2[i,totalpar] <- 0
}
  ## second group
  for(i in (nrow(data[which(mmlrval@Data$group==0),])+1):nrow(data)){
    normconst <- integrate(jointdist, -6, 6, x= as.numeric(data[i,1:ncol(data)]), 
                           alph=avec, bet=bvec, 
                           mu=personpars[1], std=sqrt(personpars[2]))$value

    for(j in 1:ncol(data)){
      par1loc <- (j-1)*2 + 1

      score2[i,par1loc] <- (1/normconst) * integrate(fun1, -6, 6, x=data[i,j],
                                                     alph=itempars[par1loc],
                                                     bet=itempars[(par1loc+1)],
                                                     xvec = 
                                                     as.numeric(data[i,1:ncol(data)]), 
                                                     avec=avec,
                                                     bvec=bvec, 
                                                     mu=personpars[1], 
                                                     std=
                                                     sqrt(personpars[2]))$value

      score2[i,(par1loc+1)] <- (1/normconst) * integrate(fun2, -6, 6, 
                                                         x=data[i,j],
                                                         alph=itempars[par1loc],
                                                         bet=
                                                         itempars[(par1loc+1)],
                                                         xvec = 
                                                         as.numeric(data[i,1:ncol(data)]), 
                                                         avec=avec, bvec=bvec, 
                                                         mu=personpars[1], 
                                                         std=sqrt(personpars[2]))$value
     }
    score2[i,(totalpar-1)] <- (1/normconst) * integrate(fun3, -6, 6, x=data[i,j],
                                             alph=itempars[par1loc],
                                             bet=itempars[(par1loc+1)],
                                             xvec = as.numeric(data[i,1:ncol(data)]), 
                                             avec=avec,
                                             bvec=bvec, mu=personpars[1], 
                                             std=sqrt(personpars[2]))$value
    score2[i,totalpar] <- (1/normconst) * integrate(fun4, -6, 6, x=data[i,j],
                                             alph=itempars[par1loc],
                                             bet=itempars[(par1loc+1)],
                                             xvec = as.numeric(data[i,1:ncol(data)]), 
                                             avec=avec,
                                             bvec=bvec, mu=personpars[1], 
                                             std=sqrt(personpars[2]))$value
  }
  ## permutate score2 into the same parameter order.
  score <- matrix(NA, nrow(data), 2*ncol(data))
  for (j in 1:ncol(data)){
    score[,j] <- score2[,(2*j-1)]
    score[,(ncol(data)+j)] <- score2[,(2*j)]
  }
  
  score <- cbind(score, score2[,(2*ncol(data)+1):totalpar])

  return(score)
}

## only works for item parameters.
info_full.mirt <- function(mmlrval,...){
    mmlrval@information[1:(2*ncol(mmlrval@Data$data)), 1:(2*ncol(mmlrval@Data$data))]
}
