#===========================================#
#===Install and load the following packages
#===========================================#
library(gamlss)
library(evd)
library(dplyr)
library(Rlab) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#This file requires the prior execution of the GGUMBEL and RGGUMBEL files.
amostras = seq(100,1000,50)

for (nss in 1:length(amostras)) {
  set.seed(100) 
  n  <- amostras[nss]
  z1 <- rlogis(n = n) 
  z2 <- rbern(n = n,0.6)
  #========================================
  Z1 <- model.matrix(~z1+z2)
  Z2 <- as.matrix(rep(1,n),byrow=TRUE,ncol=1)
  Z3 <- as.matrix(rep(1,n),byrow=TRUE,ncol=1)
  beta1     <- c(2.183,0.013,0.035)
  beta2     <- c(3.521) 
  beta3     <- c(1.305) 
  mu        <- c(exp(Z1%*%beta1))
  sigma     <- c(Z2%*%beta2)
  nu        <- c(exp(Z3%*%beta3))
  estimados <- c()
  see       <- c()
  #=========================================
  replicas  <- 3000 
  for(i in 1 :replicas){
    flag <- 0
    while(flag==0){
      tryCatch({
        tau   <- rep(0.1,n)
        y     <- rRGGUMBEL(n, mu = mu, sigma = sigma, nu = nu, tau = tau)
        aux   <- gamlss(y~1, sigma.formula=~1, nu.formula=~1, family = GGUMBEL, n.cyc = 2000)
        aux7  <- gamlss(y~z1+z2, sigma.formula=~1, nu.formula=~1, family = RGGUMBEL, mu.start = aux$mu.fv, sigma.start = aux$sigma.fv, nu.start = aux$nu.fv, tau.start = 0.5, tau.fix = TRUE, n.cyc = 2000)  
        aux6  <- gamlss(y~z1+z2, sigma.formula=~1, nu.formula=~1, family = RGGUMBEL, mu.start = aux7$mu.fv, sigma.start = aux7$sigma.fv, nu.start = aux7$nu.fv, tau.start = 0.4, tau.fix = TRUE, n.cyc = 2000)  
        aux5  <- gamlss(y~z1+z2, sigma.formula=~1, nu.formula=~1, family = RGGUMBEL, mu.start = aux6$mu.fv, sigma.start = aux6$sigma.fv, nu.start = aux6$nu.fv, tau.start = 0.3, tau.fix = TRUE, n.cyc = 2000)  
        aux4  <- gamlss(y~z1+z2, sigma.formula=~1, nu.formula=~1, family = RGGUMBEL, mu.start = aux5$mu.fv, sigma.start = aux5$sigma.fv, nu.start = aux5$nu.fv, tau.start = 0.2, tau.fix = TRUE, n.cyc = 2000)  
        aux3  <- gamlss(y~z1+z2, sigma.formula=~1, nu.formula=~1, family = RGGUMBEL, mu.start = aux4$mu.fv, sigma.start = aux4$sigma.fv, nu.start = aux4$nu.fv, tau.start = 0.1, tau.fix = TRUE, n.cyc = 2000)  
        logL  <- gen.likelihood(aux3) 
        logL()
        logL(c(coef(aux3,"mu"), coef(aux3,"sigma"),coef(aux3,"nu"), tau = 0.1))
        H   <- optimHess(c(coef(aux3,"mu"),coef(aux3,"sigma"),coef(aux3,"nu"),tau = 0.1), logL)
        se  <- sqrt(diag(solve(H[1:5,1:5])))
        summary(aux3)
        if((aux3$converged==1)&(sum(se=="NA")==0)&(sum(se=="NaN")==0)&(sum(se=="Error")==0))
        {
          estimados <- rbind(estimados,c(aux3$mu.coefficients,aux3$sigma.coefficients,aux3$nu.coefficients))			
          see       <- rbind(see,se)
          flag<-1
          print(i)
        }else{flag<-0}  
        
      },error = function(e) flag = 0) #trycatch
    } #while
  }#for
  write.table(estimados,paste("1ESTIMADOS",n,".txt"), sep="\t",row.names=FALSE)
  write.table(see,paste("1SE",n,".txt"),sep="\t",row.names=FALSE)
}
