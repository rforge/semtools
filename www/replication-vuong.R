####################
## Preliminaries  ##
####################

## packages
library("lavaan", quietly=TRUE)
library("mvtnorm", quietly=TRUE)
library("lattice", quietly=TRUE)
library("dr", quietly=TRUE)
library("ggplot2")
library("grid") # for unit
library("mnormt")
library("data.table")
library("simsem")
library("boot")
library("tables")

## auxiliary code
source("vuong_source.R")
source("sim-vuong.R")

## real-world data
load("burnout.rda")


#################
## Application ##
#################

## model definitions
model4 <- "
    # regression
    PA ~ RA + DP + EE + SE + ELC
    DP ~ RC + EE + CC
    EE ~ RC + CC + SE
    SE ~ DM + SS + PS
    ELC ~ DM + SE
    # latent variables
    RA =~ RA1 + RA2
    RC =~ RC1 + RC2 + WO1 + WO2
    CC =~ CC1 + CC2 + CC3 + CC4
    DM =~ DM1 + DM2
    SS =~ SS1 + SS2
    PS =~ PS1 + PS2
    SE =~ SE1 + SE2 + SE3
    ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
    EE =~ EE1 + EE2 + EE3
    DP =~ DP1 + DP2
    PA =~ PA1 + PA2 + PA3
"
model4_fit <- cfa(model4, data=burnout, meanstructure=TRUE)

model5 <- "
    # regression
    PA ~ RA + DP + EE + SE + ELC
    DP ~ RC + EE + CC + SE
    EE ~ RC + CC + SE
    SE ~ DM + SS + PS
    ELC ~ DM + SE
    # latent variables
    RA =~ RA1 + RA2
    RC =~ RC1 + RC2 + WO1 + WO2
    CC =~ CC1 + CC2 + CC3 + CC4
    DM =~ DM1 + DM2
    SS =~ SS1 + SS2
    PS =~ PS1 + PS2
    SE =~ SE1 + SE2 + SE3
    ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
    EE =~ EE1 + EE2 + EE3
    DP =~ DP1 + DP2
    PA =~ PA1 + PA2 + PA3
"
model5_fit <- cfa(model5, data=burnout, meanstructure=TRUE)

model7_cfa <- "
  # regression
    PA ~ RA + DP + SE
    DP ~ EE + CC + SE
    EE ~ RC + CC + SE
    SE ~ DM + SS
    ELC ~ SE + RC
  # latent variables
    RA =~ RA1 + RA2
    RC =~ RC1 + RC2 + WO1 + WO2
    CC =~ CC1 + CC2 + CC3 + CC4
    DM =~ DM1 + DM2
    SS =~ SS1 + SS2
    PS =~ PS1 + PS2
    SE =~ SE1 + SE2 + SE3
    ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
    EE =~ EE1 + EE2 + EE3
    DP =~ DP1 + DP2
    PA =~ PA1 + PA2 + PA3
  # covariance
    ELC ~~ 0*PA
"
model7_fit <- cfa(model7_cfa, data=burnout, meanstructure=TRUE)

## NET
# model5 vs. model7 
cov7 <- unclass(fitted(model7_fit)$cov)
mean7 <- unclass(fitted(model7_fit)$mean)
model5_cov7 <- cfa(model5, sample.cov=cov7, sample.mean=mean7, sample.nobs=599, meanstructure=TRUE)
net57 <- model5_cov7 @ Fit @ test[[1]][["stat"]] # model7 is not nested in model5

# model4 vs. model7
model4_cov7 <- cfa(model4, sample.cov=cov7, sample.mean=mean7, sample.nobs=599, meanstructure=TRUE)
net47 <- model4_cov7 @ Fit @ test[[1]][["stat"]] # model7 is not nested in model4

## Vuong tests
burn47 <- vuong(model4_fit, model7_fit, burnout)
burn57 <- vuong(model5_fit, model7_fit, burnout)


##################
## Simulation 1 ##
##################

if(file.exists("sim1.rda")){
  load("sim1.rda")
} else {
  sim1 <- simulation(parval = seq(0, .5, .1))
  levels(sim1$test) <- c("Overlap","LRT","CondLRT","Pref","BIC","varBIC","NET","NETval")
  save(sim1, file="sim1.rda")
}

## Figure 3
g <- ggplot(subset(sim1, test %in% c("Overlap","LRT","BIC","NET")),
            aes(x=parval, y=power, group=test))
g + geom_point(aes(shape=test), size=1.8, alpha=0.8) +
    geom_line(aes(linetype=test), size=.45, alpha=0.8) +
    facet_grid(. ~ nobs, labeller=label_bquote(n==.(x))) + theme_bw() +
    xlab("d") + ylab(expression(paste("P(", M[A] ," favored)"))) +
    theme(axis.title.x=element_text(size=rel(0.8), face='italic'),
          axis.title.y=element_text(size=rel(0.8)),
          axis.text=element_text(size=rel(0.6)),
          axis.ticks.length=unit(0,"in"),
          legend.key=element_rect(colour="white"),
          legend.text=element_text(size=rel(0.6)),
          legend.title=element_blank(),
          panel.margin=unit(0.2,"lines"),
          strip.text.x=element_text(size=rel(0.8)),
          strip.background=element_rect(color="white",fill="white"))


##################
## Simulation 2 ##
##################

if(file.exists("sim2.rda")){
  load("sim2.rda")
} else {
  sim2 <- simulation2()
  save(sim2, file="sim2.rda")
}


## Table 1
## select Vuong and boot results
dd <- sim2[grepl("Vuong|boot", sim2$test), ]

## create a factor variable of method (Vuong, boot)
dd$method <- factor(grepl("Vuong", dd$test),
                    levels=c(TRUE, FALSE), labels=c("Vuong", "Boot"))

## create a factor variable of models
dd$model <- factor(ifelse(grepl("AB", dd$test), "A-B",
                          ifelse(grepl("BC", dd$test), "B-C", "C-A")))

## delete "_" to suppress warning on "$" compiling latex table
dd$test <- factor(gsub("_", " ", gsub("_Vuong|_boot|_AB|_BC|_CA", "", dd$test)), labels=c("Avg Width","Endpoint SD","Coverage"))

## using 'tables' package
Tab <- tabular(Justify(c)*Heading(Models)*model * Justify(c)*Heading("$n$")*nobs ~
               Format(digits=3)*
               Heading()*test * Heading()*method *
               Heading()*(result) * Heading()*(identity),
               data = dd)

## apply booktab style to all following table with 'latex' function in Hmisc
booktab_mode <- booktabs()
latex(Tab)
