simordl2bb <- function(freq, nproc, nrep, ...){
  ## assumes freq includes the final cumulative
  ## proportion of 1:
  nbin <- length(freq)-1
  bp <- cumsum(freq)[1:nbin]
  res <- matrix(NA, nrep, nbin)

  muvec <- rep(0, nbin)
  ## find multivariate normal covariance matrix between bins
  cmat <- matrix(NA, nbin, nbin)
  for(i in 1:nbin){
    for(j in i:nbin){
      tmp <- sqrt(bp[i]*(1-bp[j]))/sqrt((1-bp[i])*bp[j])
      cmat[j,i] <- tmp
      cmat[i,j] <- tmp
    }
  }

  ## create block diagonal matrix so we can simulate all
  ## parameters at once
  fullcmat <- matrix(0, nbin*nproc, nbin*nproc)
  for(j in 1:nproc){
    fullcmat[(nbin*(j-1) + 1):(nbin*j),(nbin*(j-1) + 1):(nbin*j)] <- cmat
  }

  draws <- rmvnorm(nrep, rep(0,nbin*nproc), fullcmat)

  ## now get L2 norm sequentially for each bin:
  ## columns 1,(nbin+1),(2nbin+1),...
  ## columns 2,(nbin+2),(2nbin+2),...
  ## TODO if we also loops through values of colnums,
  ##      (first=1, then=2, ..., nproc)
  ##      we can simulate multiple values of nproc
  colnums <- 1:nproc
  ssbin <- matrix(NA, nrep, nbin)
  for(j in 1:nbin){
    ssbin[,j] <- apply(as.matrix(draws[,(colnums-1)*nbin + j]), 1,
                       function(x) sum(x^2))
  }

  ## now get max across the above matrix
  maxval <- apply(ssbin, 1, max)
  
  ## distribution of maximum values
  maxval
}

if(FALSE){
  library(strucchange)
  ## Load critical values that were simulated from
  ## the psychometrika paper:
  load("~/Projects/semtools/ordinal/critvals.rda")
  
  ## Define quantiles, compare to new simulation
  qs <- (1:19)/20
  tmp <- simordl2bb(rep((1/8),8), 3, 50000)
  oldsim <- cval[[2]]$computeCritval(1-qs,3)
  newsim <- quantile(tmp, qs)
  plot(oldsim, newsim)
  cbind(oldsim, newsim)

  ## 20 parameters
  tmp <- simordl2bb(rep((1/20),20), 3, 50000)
  oldsim <- cval[[4]]$computeCritval(1-qs,3)
  newsim <- quantile(tmp, qs)
  plot(oldsim, newsim)
  cbind(oldsim, newsim)

  ## Try unequal bin sizes; do the old simulation
  ## to compare (takes awhile)
  tmp <- simordl2bb(c(.2,.5,.15,.15), 3, 50000)
  tmp2 <- ordL2BB(c(.2,.5,.15,.15), nobs=5000, nproc=3, nrep=50000)
  oldsim <- tmp2$computeCritval(1-qs,3)
  newsim <- quantile(tmp, qs)
  plot(oldsim, newsim)
  cbind(oldsim, newsim)

  ## should match ordwmax when nproc==1
  ## seems to be a slight difference, but closer
  ## if there are equal bin proportions
  tmp <- simordl2bb(c(.2,.5,.15,.15), 1, 50000)
  tmp2 <- ordwmax(c(.2,.5,.15,.15))
  oldsim <- tmp2$computeCritval(.05,1)
  newsim <- quantile(tmp, .95)
  cbind(oldsim^2, newsim)
}
