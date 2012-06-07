library("lavaan")
library("strucchange")
source("../www/estfun-lavaan.R")
source("efpFunctional-cat.R")

## Define 15-category grouping variable
grp <- c(rep(1:15, each=20),15)

## Fit a model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
     
fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure=TRUE)

fit.gefp <- gefp(fit, fit = NULL, order.by = grp, sandwich = FALSE,
                 parm=7:15)

## Takes about 5min on my machine, with this warning:
#### $`data length [6772] is not a sub-multiple or multiple of the number of rows [753]`
#### matrix(rnorm(nproc * nobs), ncol = nproc)
critval <- ordL2BB(fit.gefp)

## No warnings with this, so I think the issue lies on
## line 36 of efpFunctional-cat.R
critval2 <- ordL2BB(fit.gefp, nobs = 301)
