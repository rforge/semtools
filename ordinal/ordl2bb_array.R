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
  rval <- vector(mode = "list", length = nbin)
  for(j in 1L:nbin){
    cols <- (colnums - 1L) * nbin + j
    res <- t(apply(as.matrix(draws[, cols])^2, 1, cumsum))
    if(maxproc==1) res <- matrix(res, nrep, 1)
    rval[[j]] <- res[, nproc]
  }
  rval <- do.call("pmax", rval)
  rval
}

if(FALSE){
  library("strucchange")
  library("mvtnorm")
  source("ordl2bb_array.R")
  ## Load critical values that were simulated from
  ## the psychometrika paper:
  load("~/Projects/semtools/ordinal/critvals.rda")
  
  ## Define quantiles, compare to new simulation
  qs <- c((1:19)/20, .975, .99)
  system.time(tmp <- simordl2bb(rep((1/8),8), 3, 50000))
  ## Compare to looped code
  source("ordl2bb.R")
  system.time(tmp2 <- simordl2bb(rep((1/8),8), 3, 50000))
  oldsim <- cval[[2]]$computeCritval(1-qs,3)
  newsim <- quantile(tmp, qs)
  newsim2 <- quantile(tmp2, qs)
  cbind(oldsim, newsim, newsim2)

  ## 20 parameters
  source("ordl2bb_array.R")
  system.time(tmp <- simordl2bb(rep((1/20),20), 3, 50000))
  source("ordl2bb.R")
  system.time(tmp2 <- simordl2bb(rep((1/20),20), 3, 50000))
  oldsim <- cval[[4]]$computeCritval(1-qs,3)
  newsim <- quantile(tmp, qs)
  newsim2 <- quantile(tmp2, qs)
  cbind(oldsim, newsim, newsim2)

  ## should match ordwmax when nproc==1
  tmp <- simordl2bb(c(.2,.5,.15,.15), 1, 50000)
  tmp2 <- ordwmax(c(.2,.5,.15,.15))
  oldsim <- tmp2$computeCritval(.05,1)
  newsim <- quantile(tmp, .95)
  cbind(oldsim^2, newsim)

  ## multiple values of nproc: 50000 takes awhile
  ## for nproc=1:20, but could probably get away with a
  ## smaller value of nrep.
  source("ordl2bb_array.R")
  system.time(tmp <- simordl2bb(rep((1/20),20), 1:20, 50000))
  source("ordl2bb.R")
  system.time(tmp2 <- simordl2bb(rep((1/20),20), 1:20, 50000))
  ## Comparison for each value of nproc
  for(i in 1:20){
    print(rbind(quantile(tmp[,i], qs[15:19]), quantile(tmp2[,i], qs[15:19]), cval[[4]]$computeCritval(1-qs[15:19], i)))
  }

  ## Try irregular values of nproc
  npcs <- c(3,7,9)
  source("ordl2bb_array.R")
  system.time(tmp <- simordl2bb(rep((1/20),20), npcs, 100000))
  source("ordl2bb.R")
  system.time(tmp2 <- simordl2bb(rep((1/20),20), npcs, 100000))
  for(i in 1:length(npcs)){
    print(rbind(quantile(tmp[,i], qs[15:19]), quantile(tmp2[,i], qs[15:19]), cval[[4]]$computeCritval(1-qs[15:19], npcs[i])))
  }
}
