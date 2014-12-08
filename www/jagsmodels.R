## Upon sourcing, this file creates jags syntax files for models from the paper

###########################################################
## 2f.jag
cat('
## JAGS code for preference to sure thing/risky thing. 
## 2 factor model.
## model specification 
model{
  for (i in 1:N){
    factorstar[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-int[j]+factorstar[i,1]*lambda1star[j] +
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
    Inv_sig_e[j]~dgamma(0.01,0.01) 
    int[j]~dnorm(0,0.001)
  }
   
## allow factors correlation.
  Inv_cov[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0

## transform the parameters: lambda, error variance, 
## phi(factor correlation)
  for (j in 1:J){
    Sig_e[j]<-1/Inv_sig_e[j]
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
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-lambdastar[j]*factorstar[i]+int[j]
    }
    factorstar[i] ~ dnorm(0,Inv_sig_f1) 
  }

   #Prior distribution
   for (j in 1:J){
     int[j] ~ dnorm(0,0.001)
     lambdastar[j] ~ dnorm(0,1)
   }

   for (j in 1:J){
     Inv_sig_e[j] ~ dgamma(0.01,0.01)
   }

   Inv_sig_f1 ~ dgamma(0.01,0.01)

   #Transform the variance parameters
   for (j in 1:J){
     Sig_e[j]<-1/Inv_sig_e[j]
     lambda[j]<-ifelse(lambdastar[1]<0,-1,1)*
                  lambdastar[j]*sqrt(1/Inv_sig_f1)
   }

}', file="1f.jag")

###########################################################
## 2fre.jag
cat('
## JAGS code for preference to sure thing/risky thing. 2factor 
## reduced model. orthogonal in Inv_cov. 

#model specification 
model{
  for (i in 1:N){
    factorstar[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
     for (j in 1:J){
       y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
       MuY[i,j]<-int[j]+factorstar[i,1]*lambda1star[j]+
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
    Inv_sig_e[j]~dgamma(0.01,0.01) 
    int[j]~dnorm(0,0.001)
  }
  
  #Allow factors correlation.
  Inv_cov[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0

  #Transform the parameters
  for (j in 1:J){
    Sig_e[j]<-1/Inv_sig_e[j]
    lambda1[j]<-ifelse(lambda1star[1]<0,-1,1)*lambda1star[j]*
                  sqrt(Cov[1,1])
    lambda2[j]<-ifelse(lambda2star[6]<0,-1,1)*lambda2star[j]*
                  sqrt(Cov[2,2])
  }
    phi <-ifelse((lambda1star[1]*lambda2star[6])<0,-1,1)*
            Cov[1,2]/sqrt(Cov[1,1]*Cov[2,2])

}', file="2fre.jag")

###########################################################
## mi_loading.jag
cat('
## mi_loading example; change filename.jag on the last line

model{
## Frame==TRUE 
  for (i in 1:(N/2-1)){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2])
    Mu[i,1]<-0
    Mu[i,2]<-0  
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-int[j]+factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j]
    }
  }
## Frame==FALSE
  for (i in N/2:N){
    factor2[i,1:2]~dmnorm(Mu2[i,1:2], Inv_cov2[1:2,1:2])
    Mu2[i,1]<-0
    Mu2[i,2]<-0  
    for (j in 1:J){
      y[i,j]~dnorm(MuY2[i,j], Inv_sig_e2[j])
      MuY2[i,j]<-int2[j]+factor2[i,1]*(lambda1[j]+t*delta1[j]) + 
                 factor2[i,2]*(lambda2[j]+t*delta2[j])
    }
  }

## prior distribution, can be changed to use informative prior
## prior for lambda1 and delta1
  lambda1[1]<-1
  delta1[1]<-0

  for (j in 2:J/2){
    lambda1[j]~dnorm(0,0.001)
    delta1[j]~dnorm(0,0.001)
  }

  for (j in (J/2+1):J){
    lambda1[j]<-0
    delta1[j]<-0
  }
## prior for lambda2 and delta2
  for (j in 1:J/2){
    lambda2[j]<-0
    delta2[j]<-0
  }

  lambda2[J/2+1]<-1
  delta2[J/2+1]<-0
 
  for (j in (J/2+2):J){
    lambda2[j]~dnorm(0,0.001)
    delta2[j]~dnorm(0,0.001)
   }

## prior for Inv_sig_e, Inv_sig_ee, int and int2
  for (j in 1:J){
    Inv_sig_e[j]~dgamma(0.01,0.01)
    Inv_sig_e2[j]~dgamma(0.01,0.01)
    int[j]~dnorm(0,0.001) 
    int2[j]~dnorm(0,0.001) 
  }
   
## allow factors correlation.
  Inv_cov[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0
  Inv_cov2[1:2,1:2]~dwish(r[1:2,1:2],2)
  Cov2[1:2,1:2]<-inverse(Inv_cov2[1:2, 1:2])
  r[1,1]<-1
  r[2,2]<-1
  r[1,2]<-0
  r[2,1]<-0


}', file="mi_loading.jag")

###########################################################
## mi_loadingintercept.jag
cat('
## JAGS code for preference to sure thing/risky thing. Link model:
## measurement invariance for loading+intercept w.r.t. Frame. 

## model specification 
model{
## Frame==TRUE 
  for (i in 1:(N/2-1)){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2])
    Mu[i,1]<-0
    Mu[i,2]<-0  
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-int[j]+factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j]
    }
  }
## Frame==FALSE
  for (i in N/2:N){
    factor2[i,1:2]~dmnorm(Mu2[i,1:2], Inv_cov2[1:2,1:2])
    Mu2[i,1]<-0
    Mu2[i,2]<-0
    for (j in 1:J){
      y[i,j]~dnorm(MuY2[i,j], Inv_sig_e2[j])
      MuY2[i,j]<- (int[j] + t*delta[j])+
                  factor2[i,1]*(lambda1[j]+t*delta1[j]) +
                  factor2[i,2]*(lambda2[j]+t*delta2[j])
    }
  }

## prior distribution, can be changed to use informative prior
## prior for lambda1 and delta1
  lambda1[1]<-1
  delta1[1]<-0

  for (j in 2:J/2){
    lambda1[j]~dnorm(0,0.001)T(0,)
    delta1[j]~dnorm(0,0.001)T(0,)
  }

  for (j in (J/2+1):J){
    lambda1[j]<-0
    delta1[j]<-0
  }
## prior for lambda2 and delta2
  for (j in 1:J/2){
    lambda2[j]<-0
    delta2[j]<-0
  }

  lambda2[J/2+1]<-1
  delta2[J/2+1]<-0
 
  for (j in (J/2+2):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
    delta2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e, Inv_sig_e2, int2 and delta
  for (j in 1:J){
    Inv_sig_e[j]~dgamma(0.01,0.01)
    Inv_sig_e2[j]~dgamma(0.01,0.01)
    int[j]~dnorm(0,0.001) 
    delta[j]~dnorm(0,0.001) 
  }


   
## allow factors correlation.
  Inv_cov[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0
  Inv_cov2[1:2,1:2]~dwish(r[1:2,1:2],2)
  Cov2[1:2,1:2]<-inverse(Inv_cov2[1:2, 1:2])
  r[1,1]<-1
  r[2,2]<-1
  r[1,2]<-0
  r[2,1]<-0

}', file="mi_loadingintercept.jag")

###########################################################
## 2freg.jag
cat('
## JAGS code for the impact of Frame, Num and Frame Num 
## interaction on the mean of 2 factor model

## model specification 
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:(J/2)){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + beta1[1]*Frame[i] +
                beta2*Num[i] + beta3*Num[i]*Frame[i]
    }
    for (j in (J/2+1):J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + beta1[2]*Frame[i] +
                beta2*Num[i] + beta3*Num[i]*Frame[i]
    }
  }
## prior distribution, can be changed to use informative prior

## prior for lambda1
  for (j in 1:(J/2)){
    lambda1[j]~dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j]<-0
  }

## prior for lambda2
  for (j in 1:(J/2)){
    lambda2[j]<-0
  }
  for (j in (J/2+1):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, beta2 and beta3
   for (j in 1:J){
     Inv_sig_e[j]~dgamma(0.01,0.01) 
     beta0[j]~dnorm(0,0.001)
   }
   beta1[1]~dnorm(0,0.001)
   beta1[2]~dnorm(0,0.001)
   beta2~dnorm(0,0.001)
   beta3~dnorm(0,0.001)

## allow factors correlation.
  Inv_cov[1:2,1:2]<-inverse(Cov[1:2, 1:2])
  Cov[1,1]<-1
  Cov[2,2]<-1
  Cov[1,2]<-phi
  Cov[2,1]<-phi

## prior for correlation between factors
  phi~dunif(-1,1)

}', file="2freg.jag")

###########################################################
## linkinteraction.jag
cat('
## Link model for Num*Frame interaction term.  

## model specification 
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:(J/2)){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + beta1[1]*Frame[i] +
                beta2*Num[i] + t*beta3*Num[i]*Frame[i]
    }
    for (j in (J/2+1):J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + beta1[2]*Frame[i] +
                beta2*Num[i]
    }
  }
## prior distribution, can be changed to use informative prior

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j]~dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j]<-0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j]<-0
  }
  for (j in (J/2+1):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, beta2 and beta3
   for (j in 1:J){
     Inv_sig_e[j]~dgamma(0.01,0.01) 
     beta0[j]~dnorm(0,0.001)
   }
   beta1[1]~dnorm(0,0.001)
   beta1[2]~dnorm(0,0.001)
   beta2~dnorm(0,0.001)
   beta3~dnorm(0,0.001)

## allow factors correlation.
  Inv_cov[1:2,1:2]<-inverse(Cov[1:2, 1:2])
  Cov[1,1]<-1
  Cov[2,2]<-1
  Cov[1,2]<-phi
  Cov[2,1]<-phi

## prior for correlation between factors
  phi~dunif(-1,1)
    
}', file="linkinteraction.jag")

###########################################################
## linkFramesure1.jag
cat('
## Frame term on attractiveness, one parameter. 

## model specification 
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:J/2){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + t*beta1*Frame[i] + 
                beta2*Num[i]
    }
    for (j in (J/2+1):J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + beta2*Num[i]
    }
  }
## prior distribution, can be changed to use informative prior

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j]~dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j]<-0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j]<-0
  }
  for (j in (J/2+1):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, beta2
   for (j in 1:J){
     Inv_sig_e[j]~dgamma(0.01,0.01) 
     beta0[j]~dnorm(0,0.001)
   }
   beta1~dnorm(0,0.001)
   beta2~dnorm(0,0.001)

## allow factors correlation.
  Inv_cov[1:2,1:2]<-inverse(Cov[1:2, 1:2])
  Cov[1,1]<-1
  Cov[2,2]<-1
  Cov[1,2]<-phi
  Cov[2,1]<-phi

## prior for correlation between factors
  phi~dunif(-1,1)

}', file="linkFramesure1.jag")

###########################################################
## linkFramerisky1.jag
cat('
## Frame term on attractiveness, one parameter. 

## model specification 
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0
    for (j in 1:J/2){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + 
                beta2*Num[i]
    }
    for (j in (J/2+1):J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-beta0[j]+ factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j] + t*beta1*Frame[i] + 
                beta2*Num[i]
    }
  }
## prior distribution, can be changed to use informative prior

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j]~dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j]<-0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j]<-0
  }
  for (j in (J/2+1):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, beta2
   for (j in 1:J){
     Inv_sig_e[j]~dgamma(0.01,0.01) 
     beta0[j]~dnorm(0,0.001)
   }
   beta1~dnorm(0,0.001)
   beta2~dnorm(0,0.001)

## allow factors correlation.
  Inv_cov[1:2,1:2]<-inverse(Cov[1:2, 1:2])
  Cov[1,1]<-1
  Cov[2,2]<-1
  Cov[1,2]<-phi
  Cov[2,1]<-phi

## prior for correlation between factors
  phi~dunif(-1,1)
}', file="linkFramerisky1.jag")

###########################################################
## 2fpred.jag
cat('
## JAGS code for the impact of Frame, Num and Frame Num interaction on 
## the mean of 2 factor model and on the actual choice

## model specification
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0

## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:(J/2)){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] +
                  factor[i,2]*lambda2[j] + beta1[1]*Frame[i] + 
                  beta2*Num[i] + beta3*Frame[i]*Num[i]
    }
    for (j in (J/2+1):J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] +
                  factor[i,2]*lambda2[j] + beta1[2]*Frame[i] + 
                  beta2*Num[i] + beta3*Frame[i]*Num[i]
    }
## impact of Frame, Num and Frame Num interaction on choice
    for (k in 1:K){  
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + b1*Frame[i] + 
                  b2*(Num[i]-Numc) + b3*Frame[i]*(Num[i]-Numc) + 
                  g1*MuY[i,k+J/2]+ g2*MuY[i,k] + 
                  g3*MuY[i,k+J/2]*(Num[i]-Numc) + 
                  g4*MuY[i,k]*(Num[i]-Numc) 
    }
}

## prior distribution, can be changed to use informative prior
## prior for, beta2~beta3, b1-b3, g1-g4
  beta2 ~ dnorm(0,0.001)
  beta3 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)
  b3 ~ dnorm(0,0.001)
  g1 ~ dnorm(0,0.001)
  g2 ~ dnorm(0,0.001)
  g3 ~ dnorm(0,0.001)
  g4 ~ dnorm(0,0.001)

## prior for b0 and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    Inv_sig_ee[k] ~ dgamma(0.01,0.01)
  }

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda2[j] ~ dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, b1
  for (j in 1:J){
    Inv_sig_e[j] ~ dgamma(0.01,0.01) 
    beta0[j] ~ dnorm(0,0.001)
  }
  beta1[1] ~ dnorm(0,0.001)
  beta1[2] ~ dnorm(0,0.001)

## allow factors correlation.
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi

## prior for correlation between factors
  phi ~ dunif(-1,1)
}', file="2fpred.jag")

###########################################################
## linkinteractionchoicesure.jag
cat('
## Link model for sure*Num interaction term on choice model.

## model specification
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0

## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] + 
                  factor[i,2]*lambda2[j] 
    }
## impact of Frame, Num and Frame Num interaction on choice
    for (k in 1:K){  
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + b1*Frame[i] + b2*(Num[i]-Numc) +
                  g1*MuY[i,k+J/2]+ g2*MuY[i,k] + 
                  g3*MuY[i,k+J/2]*(Num[i]-Numc) + 
                  t*g4*MuY[i,k]*(Num[i]-Numc) 
    }
}

## prior distribution, can be changed to use informative prior
## prior for beta2-beta3, b1-b3, g1-g4
  beta2 ~ dnorm(0,0.001)
  beta3 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)
  b3 ~ dnorm(0,0.001)
  g1 ~ dnorm(0,0.001)
  g2 ~ dnorm(0,0.001)
  g3 ~ dnorm(0,0.001)
  g4 ~ dnorm(0,0.001)

 
## prior for b0 and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    Inv_sig_ee[k] ~ dgamma(0.01,0.01)
  }


## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda2[j] ~ dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, b1
  for (j in 1:J){
    Inv_sig_e[j] ~ dgamma(0.01,0.01) 
    beta0[j] ~ dnorm(0,0.001)
    beta1[j] ~ dnorm(0,0.001)
  }

## allow factors correlation.
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi

## prior for correlation between factors
  phi ~ dunif(-1,1)
}', file="linkinteractionchoicesure.jag")

###########################################################
## linkinteractionchoicerisky.jag
cat('
## Link model for risky*Num interaction term on choice model.

## model specification
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0

## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] + 
                  factor[i,2]*lambda2[j] 
    }
## impact of Frame, Num and Frame Num interaction on choice
    for (k in 1:K){  
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + b1*Frame[i] + b2*(Num[i]-Numc) +
                  g1*MuY[i,k+J/2]+ g2*MuY[i,k] + 
                  t*g3*MuY[i,k+J/2]*(Num[i]-Numc) + 
                  g4*MuY[i,k]*(Num[i]-Numc) 
    }
}

## prior distribution, can be changed to use informative prior
## prior for b1-b3, g1-g4.
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)
  b3 ~ dnorm(0,0.001)
  g1 ~ dnorm(0,0.001)
  g2 ~ dnorm(0,0.001)
  g3 ~ dnorm(0,0.001)
  g4 ~ dnorm(0,0.001)

 
## prior for b0 and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    Inv_sig_ee[k] ~ dgamma(0.01,0.01)
  }

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda2[j] ~ dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, b1
  for (j in 1:J){
    Inv_sig_e[j] ~ dgamma(0.01,0.01) 
    beta0[j] ~ dnorm(0,0.001)
    beta1[j] ~ dnorm(0,0.001)
  }

## allow factors correlation.
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi

## prior for correlation between factors
  phi ~ dunif(-1,1)
}', file="linkinteractionchoicerisky.jag")

###########################################################
## linkFramechoice1.jag
cat('
## Link model for Frame term on choice.

## model specification
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0

## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] + 
                  factor[i,2]*lambda2[j] 
    }
## impact of Frame, Num and Frame Num interaction on choice
    for (k in 1:K){  
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + t*b1*Frame[i] + b2*(Num[i]-Numc) +
                  g1*MuY[i,k+J/2]+ g2*MuY[i,k] + 
                  g3*MuY[i,k+J/2]*(Num[i]-Numc) + 
                  g4*MuY[i,k]*(Num[i]-Numc) 
    }
}

## prior distribution, can be changed to use informative prior
## prior for b1-b2, g1-g4
  b2 ~ dnorm(0,0.001)
  g1 ~ dnorm(0,0.001)
  g2 ~ dnorm(0,0.001)
  g3 ~ dnorm(0,0.001)
  g4 ~ dnorm(0,0.001)

 
## prior for b0 and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    Inv_sig_ee[k] ~ dgamma(0.01,0.01)
  }
  b1 ~ dnorm(0,0.001)

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda2[j] ~ dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, b1
  for (j in 1:J){
    Inv_sig_e[j] ~ dgamma(0.01,0.01) 
    beta0[j] ~ dnorm(0,0.001)
  }

## allow factors correlation.
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi
 
## prior for correlation between factors
  phi ~ dunif(-1,1)
}', file="linkFramechoice1.jag")

###########################################################
## linkb3choice.jag
cat('
## link model for interaction term on the choice.

## model specification
model{
  for (i in 1:N){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2]) 
    Mu[i,1]<-0
    Mu[i,2]<-0

## impact of Frame, Num and Frame Num interaction on factor mean
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j] <- beta0[j] + factor[i,1]*lambda1[j] + 
                  factor[i,2]*lambda2[j] 
    }
## impact of Frame, Num and Frame Num interaction on choice
    for (k in 1:K){  
      X[i,k] ~ dnorm(MuX[i,k], Inv_sig_ee[k])
      MuX[i,k] <- b0[k] + b1*Frame[i] + 
                  b2*(Num[i]-Numc) + 
                  t*b3*Frame[i]*(Num[i]-Numc) + 
                  g1*MuY[i,k+J/2]+g2*MuY[i,k] + 
                  g3*MuY[i,k+J/2]*(Num[i]-Numc) + 
                  g4*MuY[i,k]*(Num[i]-Numc)
    }
}

## prior distribution, can be changed to use informative prior
## prior for beta2-beta3, b1-b3, g1-g4
  beta2 ~ dnorm(0,0.001)
  beta3 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)
  b3 ~ dnorm(0,0.001)
  g1 ~ dnorm(0,0.001)
  g2 ~ dnorm(0,0.001)
  g3 ~ dnorm(0,0.001)
  g4 ~ dnorm(0,0.001)

## prior for b0 and Inv_sig_ee
  for (k in 1:K){
    b0[k] ~ dnorm(0,0.001)
    Inv_sig_ee[k] ~ dgamma(0.01,0.01)
  }

## prior for lambda1
  for (j in 1:J/2){
    lambda1[j] ~ dnorm(0,0.001)T(0,)
  }
  for (j in (J/2+1):J){
    lambda1[j] <- 0
  }

## prior for lambda2
  for (j in 1:J/2){
    lambda2[j] <- 0
  }
  for (j in (J/2+1):J){
    lambda2[j] ~ dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e and beta0, beta1, b1
  for (j in 1:J){
    Inv_sig_e[j] ~ dgamma(0.01,0.01) 
    beta0[j] ~ dnorm(0,0.001)
    beta1[j] ~ dnorm(0,0.001)
  }

## allow factors correlation.
  Inv_cov[1:2,1:2] <- inverse(Cov[1:2, 1:2])
  Cov[1,1] <- 1
  Cov[2,2] <- 1
  Cov[1,2] <- phi
  Cov[2,1] <- phi

## prior for correlation between factors
  phi ~ dunif(-1,1)

}', file="linkb3choice.jag")

###########################################################
## mi_intercept.jag
cat('
## JAGS code for preference to sure thing/risky thing. Link model:
## measurement invariance for loading+intercept w.r.t. Frame. 

## model specification 
model{
## Frame==TRUE 
  for (i in 1:(N/2-1)){
    factor[i,1:2]~dmnorm(Mu[i,1:2], Inv_cov[1:2,1:2])
    Mu[i,1]<-0
    Mu[i,2]<-0  
    for (j in 1:J){
      y[i,j]~dnorm(MuY[i,j], Inv_sig_e[j])
      MuY[i,j]<-int[j]+factor[i,1]*lambda1[j] + 
                factor[i,2]*lambda2[j]
    }
  }
## Frame==FALSE
  for (i in N/2:N){
    factor2[i,1:2]~dmnorm(Mu2[i,1:2], Inv_cov2[1:2,1:2])
    Mu2[i,1]<-0
    Mu2[i,2]<-0
    for (j in 1:J){
      y[i,j]~dnorm(MuY2[i,j], Inv_sig_e2[j])
      MuY2[i,j]<- (int[j] + t*delta[j])+
                  factor2[i,1]*lambda1[j] +
                  factor2[i,2]*lambda2[j]
    }
  }

## prior distribution, can be changed to use informative prior
## prior for lambda1 and delta1
  lambda1[1]<-1
  delta1[1]<-0

  for (j in 2:J/2){
    lambda1[j]~dnorm(0,0.001)T(0,)
    delta1[j]~dnorm(0,0.001)T(0,)
  }

  for (j in (J/2+1):J){
    lambda1[j]<-0
    delta1[j]<-0
  }
## prior for lambda2 and delta2
  for (j in 1:J/2){
    lambda2[j]<-0
    delta2[j]<-0
  }

  lambda2[J/2+1]<-1
  delta2[J/2+1]<-0
 
  for (j in (J/2+2):J){
    lambda2[j]~dnorm(0,0.001)T(0,)
    delta2[j]~dnorm(0,0.001)T(0,)
  }

## prior for Inv_sig_e, Inv_sig_e2, int2 and delta
  for (j in 1:J){
    Inv_sig_e[j]~dgamma(0.01,0.01)
    Inv_sig_e2[j]~dgamma(0.01,0.01)
    int[j]~dnorm(0,0.001) 
    delta[j]~dnorm(0,0.001) 
  }


   
## allow factors correlation.
  Inv_cov[1:2,1:2]~dwish(R[1:2,1:2],2)
  Cov[1:2,1:2]<-inverse(Inv_cov[1:2, 1:2])
  R[1,1]<-1
  R[2,2]<-1
  R[1,2]<-0
  R[2,1]<-0
  Inv_cov2[1:2,1:2]~dwish(r[1:2,1:2],2)
  Cov2[1:2,1:2]<-inverse(Inv_cov2[1:2, 1:2])
  r[1,1]<-1
  r[2,2]<-1
  r[1,2]<-0
  r[2,1]<-0

}', file="mi_intercept.jag")
