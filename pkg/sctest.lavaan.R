sctest.lavaan <- function(x, order.by, parm, stat, fnl = NULL, fplot = FALSE, ...){
  ## Carry out a score-based test of a fitted lavaan model x

  require("strucchange")
  
  ## TODO Check that the fit uses MVN likelihood fn
  
  ## Couldn't gefp() have a "type" argument like efp()?
  fp <- gefp(x, fit = NULL, order.by = order.by, vcov = info_full,
             sandwich = FALSE, parm = parm)

  ## switch between stat argument names and corresponding functionals
  ## ordwmax(fp), ordL2BB(fp), maxBB, meanL2BB, supLM(0.1), rangeBB
  ## WDMo         maxLMo       DM     CvM       maxLM,      range
  ## TODO? Allow user to modify the 0.1 of supLM?
  if(is.null(fnl)){
    fnl <- switch(stat,
                  WDMo = ordwmax(fp),
                  maxLMo = ordL2BB(fp),
                  DM = maxBB,
                  CvM = meanL2BB,
                  maxLM = supLM(0.1),
                  range = rangeBB,
                  LMuo = catL2BB(fp))
  }

  ## If stat argument is messed up
  if(is.null(fnl)) stop(paste("Argument stat=",stat," is unknown",sep=""))
  
  testres <- sctest(fp, functional = fnl)

  if(fplot) plot(fp, functional = fnl, ...)

  testres
}
  
## convenience function for applying tests to real data
info_full <- function(x, ...) solve(vcov(x) * nobs(x))
