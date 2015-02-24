simordl2bb <- function(freq, nproc, nrep, ...){
  ## assumes freq includes the final cumulative proportion of 1
  nbin <- length(freq) - 1L
  bp <- cumsum(freq)[1L:nbin]
  maxproc <- max(nproc)

  ## find multivariate normal covariance matrix between bins
  cmat <- matrix(NA, nbin, nbin)
  for(i in 1L:nbin) cmat[i:nbin, i] <- cmat[i, i:nbin] <- sqrt(bp[i] * (1 - bp[i:nbin]))/sqrt((1 - bp[i]) * bp[i:nbin])

  ## create block-diagonal matrix so we can simulate all parameters at once
  fullcmat <- matrix(0, nbin * maxproc, nbin * maxproc)
  for(j in 1L:maxproc) fullcmat[(nbin * (j - 1L) + 1L):(nbin * j), (nbin * (j - 1L) + 1L):(nbin * j)] <- cmat

  draws <- rmvnorm(nrep, sigma = fullcmat)

  ## now get L2 norm sequentially for each bin:
  ## columns 1,(nbin+1),(2nbin+1),...
  ## columns 2,(nbin+2),(2nbin+2),...
  colnums <- 1L:maxproc
  ssbins <- matrix(NA, nrep, length(nproc))
  for(i in 1L:length(nproc)){
    tmpbin <- matrix(NA, nrep, nbin)
    for(j in 1L:nbin){
      cols <- ((colnums - 1L) * nbin + j)[1:nproc[i]]
      tmpbin[, j] <- rowSums(as.matrix(draws[, cols])^2)
    }
    cmd <- paste("tmpbin[, ", 1:nbin, "]", sep="", collapse=",")
    ssbins[,i] <- eval(parse(text=paste("pmax(",cmd,")")))
  }

  ssbins
}

if(FALSE){
  library("strucchange")
  library("mvtnorm")
  source("ordl2bb.R")
  ## Load critical values that were simulated from
  ## the psychometrika paper:
  load("~/Projects/semtools/ordinal/critvals.rda")
  
  ## Define quantiles, compare to new simulation
  qs <- c((1:19)/20, .975, .99)
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
  tmp <- simordl2bb(c(.2,.5,.15,.15), 1, 50000)
  tmp2 <- ordwmax(c(.2,.5,.15,.15))
  oldsim <- tmp2$computeCritval(.05,1)
  newsim <- quantile(tmp, .95)
  cbind(oldsim^2, newsim)

  ## multiple values of nproc: 50000 takes awhile
  ## for nproc=1:20, but could probably get away with a
  ## smaller value of nrep.
  tmp <- simordl2bb(rep((1/20),20), 1:20, 50000)
  ## Comparison for each value of nproc
  for(i in 1:20){
    print(rbind(quantile(tmp[,i], qs[15:19]), cval[[4]]$computeCritval(1-qs[15:19], i)))
  }

  ## Try irregular values of nproc
  npcs <- c(3,7,9)
  tmp <- simordl2bb(rep((1/20),20), npcs, 50000)
  for(i in 1:length(npcs)){
    print(rbind(quantile(tmp[,i], qs[15:19]), cval[[4]]$computeCritval(1-qs[15:19], npcs[i])))
  }
}
