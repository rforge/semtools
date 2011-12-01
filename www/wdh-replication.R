## Package and data
library("lavaan")
data("StereotypeThreat", package = "psychotools")

## Include group variable and order correspondingly
StereotypeThreat <- transform(StereotypeThreat, group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]

## Step 2: Fix loadings across groups
m2  <- c(
'ability   =~ label(rep("load_n", 4)) * numerical  + 1 * abstract + label(rep("load_v", 4)) * verbal',
'ability   ~  0 * 1',
'ability   ~~ label(c("lvar:maj_c", "lvar:min_c", "lvar:maj_t", "lvar:min_t")) * ability',

'numerical ~  label(c("mean_n:maj_c", "mean_n:min_c", "mean_n:maj_t", "mean_n:min_t")) * 1',
'abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", "mean_a:maj_t", "mean_a:min_t")) * 1',
'verbal    ~  label(c("mean_v:maj_c", "mean_v:min_c", "mean_v:maj_t", "mean_v:min_t")) * 1',

'numerical ~~ label(c("var_n:maj_c", "var_n:min_c", "var_n:maj_t", "var_n:min_t")) * numerical',
'abstract  ~~ label(c("var_a:maj_c", "var_a:min_c", "var_a:maj_t", "var_a:min_t")) * abstract',
'verbal    ~~ label(c("var_v:maj_c", "var_v:min_c", "var_v:maj_t", "var_v:min_t")) * verbal')

## Step 2a: Free numerical loading in group min_t (minority.threat)
m2a <- m2
m2a[1] <- 'ability   =~ label(c(rep("load_n", 3), "load_n:min_t")) * numerical  + 1 * abstract + label(rep("load_v", 4)) * verbal'

## Step 3: Fix variances across groups
m3 <- m2a
m3[7:9] <- c(
'numerical ~~ label(rep("var_n", 4)) * numerical',
'abstract  ~~ label(rep("var_a", 4)) * abstract',
'verbal    ~~ label(rep("var_v", 4)) * verbal')

## Step 3a: Free numerical variance in group min_t
m3a <- m3
m3a[7] <- 'numerical ~~ label(c(rep("var_n", 3), "var_n:min_t")) * numerical'

## Step 4: Fix latent variances within conditions
m4 <- m3a
m4[3] <- 'ability   ~~ label(c("lvar:maj", "lvar:min", "lvar:maj", "lvar:min")) * ability'

## Step 5: Fix certain means, free others ?Z?
m5 <- m4
m5[c(2, 4:6)] <- c(
'ability   ~  label(c(NA, "lmean:min_c", "lmean:maj_t", "lmean:min_t")) * 1 + c(0, NA, NA, NA) * 1',
'numerical ~  label(c(rep("mean_n", 3), "mean_n:min_t")) * 1',
'abstract  ~  label(rep("mean_a", 4)) * 1',
'verbal    ~  label(rep("mean_v", 4)) * 1')

## Step 5a: Free ability mean in group maj_c
m5a <- m5
m5a[5] <- 'abstract  ~  label(c("mean_a:maj_c", rep("mean_a", 3))) * 1'

## Step 5b: Free ability mean in group min_c
m5b <- m5a
m5b[5] <- 'abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", rep("mean_a", 2))) * 1'

## Step 6: Different latent mean structure ?Z?
m6 <- m5b
m6[2] <- 'ability   ~  label(c("lmean:maj", NA, "lmean:maj", NA)) * 1 + c(NA, 0, NA, 0) * 1'

## Fit all models
m <- list(m2, m2a, m3, m3a, m4, m5, m5a, m5b, m6)
names(m) <- c("2", "2a", "3", "3a", "4", "5", "5a", "5b", "6")
fm <- lapply(m, function(f)
  lavaan(f, data = StereotypeThreat, meanstructure = TRUE, group = "group")) ## , likelihood = "wishart"))

## Extract measures of fit
tab <- t(sapply(fm, fitMeasures, c("chisq", "df", "pvalue", "rmsea", "cfi")))
tab <- rbind("1" = c(0, 0, 1, 0, 1), tab)
tab <- cbind(tab,
  delta_chisq = c(NA, abs(diff(tab[, "chisq"]))),
  delta_df = c(NA, diff(tab[, "df"])))
tab <- cbind(tab, "pvalue2" = pchisq(tab[, "delta_chisq"], abs(tab[, "delta_df"]), lower.tail = FALSE))
tab <- tab[, c(2, 1, 3, 7, 6, 8, 4, 5)]
round(tab, digits = 3)

## Observation: CFI appears to be closer to this.... ?Z?
tab[, "cfi"] * 10 - 9

## Step 7: Fix all latent variances, and all numerical means ?Z?
m7 <- m6
m7[3:4] <- c(
'ability   ~~ label(rep("lvar", 4)) * ability',
'numerical ~  label(rep("mean_n", 4)) * 1')
fm7 <- lavaan(m7, data = StereotypeThreat, meanstructure = TRUE, group = "group", likelihood = "wishart")
fitMeasures(fm7, c("chisq", "df", "pvalue", "rmsea", "cfi"))
