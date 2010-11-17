# Main R script file to carry out tests of measurement invariance
# described in Merkle & Zeileis.

# Load requisite packages
library(OpenMx)
library(strucchange)

# Read data
obs.data <- dget("data.dat")
dat <- obs.data$datmat
V <- obs.data$age

# Obtain functions for score and information matrix extraction
source("fa_score_extraction.R")

# Obtain OpenMx model specification
# (Two-group model with all parameters constrained to be equal,
#  though this could also be estimated as a single-group model.)
source("def.model.R")
fa.mod <- def.model(dat, V)

# Fit model via OpenMx
mx.res <- mxRun(fa.mod)

# Organize results for strucchange:
fa.res <- list(res=mx.res, dat=dat, V=V)

# Calculate cumulative scores for first three factor loadings:
fa.cums <- gefp(fa.res, fit=NULL, scores=score.fa, vcov=inf.fa, 
                order.by=fa.res$V, sandwich=FALSE, parm=1:3)

# Plot individual cumulative score processes:
plot(fa.cums, agg=FALSE, ylim=c(-1.8,0.3),
     ylab="Cumulative scores", xlab="Age")

# Obtain results of double-max test:
sctest(fa.cums, functional="max")

# Obtain results of Cramer-von Mises test:
sctest(fa.cums, functional="meanL2")
