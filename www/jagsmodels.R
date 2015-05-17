## Upon sourcing, this file creates jags syntax files for models from the paper

###########################################################
## 2f.jag
cat('
## JAGS code for preference to sure thing/risky thing. 
## 2 factor model.
## model specification 
model{
  for (i in 1:N){
    factorstar[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov1[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e1[j])
      MuY[i,j]<-int1[j]+factorstar[i,1]*lambda1star[j] +
                  factorstar[i,2]*lambda2star[j]
    }
  }
## prior distribution, can be changed to use informative prior
## prior for lambda1star. one of lambda1star has to be fixed to 0, 
## because factors are allowed to be correlated.
  for (j in 1:J/2){
    lambda1star[j]~dnorm(0,1)
  }
  lambda1star[J/2+1]<-0
   for (j in (J/2+2):J){
    lambda1star[j]~dnorm(0,1)
  }

## prior for lambda2star 
  lambda2star[1]<-0
  for (j in 2:J){
    lambda2star[j]~dnorm(0,1)
  }

## prior for error variances and intercepts
  for (j in 1:J){
    Inv_sig_e1[j]~dgamma(0.1,0.1) 
    int1[j]~dnorm(0,0.001)
  }
   
## allow factors correlation.
  Inv_cov1[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov1[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0

## transform the parameters: lambda, error variance, 
## phi(factor correlation)
  for (j in 1:J){
    Sig_e[j]<-1/Inv_sig_e1[j]
    lambda1[j]<-ifelse(lambda1star[1]<0,-1,1)*
                  lambda1star[j]*sqrt(Cov[1,1])
    lambda2[j]<-ifelse(lambda2star[6]<0,-1,1)*
                  lambda2star[j]*sqrt(Cov[2,2])
  }
  phi<-ifelse((lambda1star[1]*lambda2star[6])<0,-1,1)*
         Cov[1,2]/sqrt(Cov[1,1]*Cov[2,2])
}', file="2f.jag")

###########################################################
## 1f.jag
cat('
## JAGS code for preference to sure thing/risky thing. 
## 1 factor model.

## model specification 
model{
  for (i in 1:N){
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e1[j])
      MuY[i,j]<-lambda1star[j]*factorstar[i]+int1[j]
    }
    factorstar[i] ~ dnorm(0,Inv_sig_f1) 
  }

   #Prior distribution
   for (j in 1:J){
     int1[j] ~ dnorm(0,0.001)
     lambda1star[j] ~ dnorm(0,1)
   }

   for (j in 1:J){
     Inv_sig_e1[j] ~ dgamma(0.1,0.1)
   }

   Inv_sig_f1 ~ dgamma(0.1,0.1)

   #Transform the variance parameters
   for (j in 1:J){
     Sig_e[j]<-1/Inv_sig_e1[j]
     lambda1[j]<-ifelse(lambda1star[1]<0,-1,1)*
                  lambda1star[j]*sqrt(1/Inv_sig_f1)
   }

}', file="1f.jag")

###########################################################
## 2fre.jag
cat('
## JAGS code for preference to sure thing/risky thing. 2factor 
## reduced model. not orthogonal in Inv_cov. 

#model specification 
model{
  for (i in 1:N){
    factorstar[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov1[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
     for (j in 1:J){
       y[i,j]~dnorm(MuY[i,j], Inv_sig_e1[j])
       MuY[i,j]<-int1[j]+factorstar[i,1]*lambda1star[j]+
                   factorstar[i,2]*lambda2star[j]
    }
  }
  #Prior distribution, can be changed to use informative prior
  for (j in 1:J/2){
    lambda1star[j]~dnorm(0,1)
  }
  for (j in (J/2+1):J){
    lambda1star[j]<-0
  }
 
  for (j in 1:J/2){
    lambda2star[j]<-0
  }
  for (j in (J/2+1):J){
    lambda2star[j]~dnorm(0,1)
  }
 
  for (j in 1:J){
    Inv_sig_e1[j]~dgamma(0.1,0.1) 
    int1[j]~dnorm(0,0.001)
  }
  
  #Allow factors correlation.
  Inv_cov1[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov1[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0

  #Transform the parameters
  for (j in 1:J){
    Sig_e[j]<-1/Inv_sig_e1[j]
    lambda1[j]<-ifelse(lambda1star[1]<0,-1,1)*lambda1star[j]*
                  sqrt(Cov[1,1])
    lambda2[j]<-ifelse(lambda2star[6]<0,-1,1)*lambda2star[j]*
                  sqrt(Cov[2,2])
  }
    phi <-ifelse((lambda1star[1]*lambda2star[6])<0,-1,1)*
            Cov[1,2]/sqrt(Cov[1,1]*Cov[2,2])

}', file="2fre.jag")

###########################################################
## 2fsemt.jag
cat('
## JAGS code for the impact of Frame, Num and Frame Num interaction on 
## the mean of 2 factor model and on the actual choice

## model specification
model{
  for (i in 1:N){
    factor[i,1:2] ~ dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1] <- tbeta1_1*beta1[1]*Frame[i] + tbeta2_1*beta2[1]*(Num[i] - Numc) +
               tbeta3_1*beta3[1]*Frame[i]*(Num[i] - Numc)
    Mu[i,2] <- tbeta1_2*beta1[2]*Frame[i] + tbeta2_2*beta2[2]*(Num[i] - Numc) +
               tbeta3_2*beta3[2]*Frame[i]*(Num[i] - Numc)

    ## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:J){
      y[i,j] ~ dnorm(MuY[i,j], Inv_sig_e1[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] +
                  factor[i,2]*lambda2[j]
    }
    ## impact of Frame, Num and Frame Num interaction on choice
    factorch[i] ~ dnorm(Much[i], Inv_var)
    Much[i] <- tnu1*nu1*factor[i,1] + tnu2*nu2*factor[i,2] +
               tbeta4*beta4*Frame[i] + tbeta5*beta5*(Num[i]-Numc) +
               tbeta6*beta6*Frame[i]*(Num[i]-Numc) +
               tnu3*nu3*factor[i,1]*(Num[i]-Numc) +
               tnu4*nu4*factor[i,2]*(Num[i]-Numc)
    for (k in 1:K){
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + lambda3[k]*factorch[i]
    }
  }
  ## Total variance of choice lv equals 1
  Inv_var <- 1 #- nu[1]^2 - nu[2]^2

  ## prior distribution, can be changed to use informative prior
  ## prior for, beta2~beta3, b1-b3, g1-g4
  nu1 ~ dnorm(0, pnu1)
  nu2 ~ dnorm(0, pnu2)
  nu3 ~ dnorm(0, pnu3)
  nu4 ~ dnorm(0, pnu4)
  for(j in 1:2){
    beta1[j] ~ dnorm(0, pbeta1)
    beta2[j] ~ dnorm(0, pbeta2)
    beta3[j] ~ dnorm(0, pbeta3)
  }
  beta4 ~ dnorm(0, pbeta4)
  beta5 ~ dnorm(0, pbeta5)
  beta6 ~ dnorm(0, pbeta6)
  
  ## prior for b0, lambda3, and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    lambda3[k] ~ dnorm(0, plambda3)T(0,)
    Inv_sig_ee[k] ~ dgamma(0.1,0.1)
  }

  ## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,plambda1)T(0,)
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
    lambda2[j] ~ dnorm(0,plambda2)T(0,)
  }

  ## prior for Inv_sig_e and beta0, b1
  for (j in 1:J){
    Inv_sig_e1[j] ~ dgamma(0.1,0.1) 
    beta0[j] ~ dnorm(0,0.001)
  }

  ## allow factors to correlate
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi

  ## prior for correlation between factors
  phi ~ dunif(-1,1)
}', file="2fsemt.jag")
