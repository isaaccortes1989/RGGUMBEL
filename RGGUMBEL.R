#===========================================#
#===Install and load the following packages
#===========================================#
library(evd)
library(gamlss)

dRGGUMBEL = function(x,mu=1,sigma=1,nu=1,tau=0.5,log=FALSE){ 
  if(any(mu<=0)) stop(paste("mu must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  h_l <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
  lfy <- log(nu)+(nu-1)*log(x)-nu*log(mu)-log(1-pgumbel(-sigma))+log(h_l)+dgumbel(((x/mu)^nu)*h_l-sigma,log=TRUE)
  if(log== FALSE)
    fy <- exp(lfy)
  else fy <- lfy
  return(fy)
}

pRGGUMBEL = function(q,mu=1,sigma=1,nu=1,tau=0.5,lower.tail=TRUE,log.p=FALSE){
  if(any(mu<=0)) stop(paste("mu must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  h_l <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
  lGy <-  log(pgumbel(((q/mu)^nu)*h_l-sigma)-pgumbel(-sigma))-log(1-pgumbel(-sigma))
  if(log.p == FALSE)
    cdf <- exp(lGy)  
  else cdf <- lGy
  if(lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1-cdf
  return(cdf)
}

rRGGUMBEL = function(n,mu=1,sigma=1,nu=1,tau=0.5){ 
  if(any(mu<=0)) stop(paste("mu must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  U   <- runif(n)
  h_l <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
  X   <- mu*((sigma+qgumbel(U-U*pgumbel(-sigma)+pgumbel(-sigma)))/h_l)^(1/nu)
  return(X)
} 

qRGGUMBEL = function(p,mu=1,sigma=1,nu=1,tau=0.5,lower.tail=TRUE,log.p=FALSE){ 
  if(any(mu<=0)) stop(paste("mu must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))    
  if(log.p == TRUE)
    p <- exp(p)  
  else p <- p
  if(lower.tail == TRUE)
    p <- p
  else p <- 1-p
  if(any(p<0)|any(p>1)) stop(paste("p must be between 0 and 1","\n",""))
  U   <- runif(n)
  h_l <- sigma+qgumbel(tau-tau*pgumbel(-sigma)+pgumbel(-sigma))
  Q   <- mu*((sigma+qgumbel(p-p*pgumbel(-sigma)+pgumbel(-sigma)))/h_l)^(1/nu)
  return(Q)
}


RGGUMBEL<-function(mu.link = "log", sigma.link = "identity",nu.link = "log",tau.link = "identity"){
  mstats <- checklink("mu.link","RGGUMBEL",substitute(mu.link),c("inverse","log","identity","own"))
  dstats <- checklink("sigma.link","RGGUMBEL",substitute(sigma.link),c("inverse","log","identity","own"))
  vstats <- checklink("nu.link","RGGUMBEL",substitute(nu.link),c("log","identity"))
  tstats <- checklink("tau.link","RGGUMBEL",substitute(tau.link),c("identity"))
  structure(list(family = c("RGGUMBEL","Reparameterized Generalized Gumbel"), 
  parameters = list(mu=TRUE,sigma=TRUE,nu=TRUE,tau=FALSE),nopar=4,type = "Continuous", 
  mu.link = as.character(substitute(mu.link)), sigma.link = as.character(substitute(sigma.link)), 
  nu.link = as.character(substitute(nu.link)), tau.link = as.character(substitute(tau.link)),
  mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, 
  nu.linkfun = vstats$linkfun, tau.linkfun = tstats$linkfun,
  mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, 
  nu.linkinv = vstats$linkinv, tau.linkinv = tstats$linkinv,
  mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, 
  nu.dr = vstats$mu.eta,tau.dr = tstats$mu.eta,
                   dldm = function(y,mu,sigma,nu,tau){ 
                   h_l  <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1  <- exp(-(y/mu)^(nu)*h_l+sigma)
                   dldm <- -(nu/mu)+nu*(y^(nu)/mu^(nu+1))*h_l*(1-k_1)
                   return(dldm)
                 },d2ldm2 = function(y,mu,sigma,nu,tau){ 
                   h_l     <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1     <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)*h_l
                   u       <- qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   l       <- (dgumbel(-sigma)*(tau-1))/dgumbel(u)                   
                   d2ldm2  <-(nu/mu^2)-nu*(nu+1)*(y^(nu)/mu^(nu+2))*h_l*(1-k_1)-nu^2*(y^(2*nu)/mu^(2*nu+2))*h_l^2*k_1
                   return(d2ldm2)
                 },dldd = function(y,mu,sigma,nu,tau){
                   h_l  <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1  <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2  <- (dgumbel(-sigma)*(tau-1))/dgumbel(qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma)))
                   dldd <- ((1+k_2)/h_l)-(dgumbel(-sigma)/(1-pgumbel(-sigma)))+(1-(y/mu)^(nu)*(1+k_2))*(1-k_1)
                   return(dldd)
                 },d2ldd2 = function(y,mu,sigma,nu,tau){
                   h_l     <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1     <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)*h_l
                   u       <- qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   l       <- (dgumbel(-sigma)*(tau-1))/dgumbel(u)                   
                   k_3    <- dgumbel(-sigma)*((exp(sigma)-1)/(1-pgumbel(-sigma)))+(dgumbel(-sigma)/(1-pgumbel(-sigma)))^2
                   k_4    <- ((tau-1)/h_l)*(dgumbel(-sigma)*(1-exp(sigma))/dgumbel(u)-(tau-1)*dgumbel(-sigma)^2*(exp(-u)-1)/dgumbel(u)^2)
                   k_5    <- -(1/h_l^2)*(1+(tau-1)*dgumbel(-sigma)/dgumbel(u))^2
                   k_6    <- -(1-(y/mu)^(nu)*(1+((tau-1)*dgumbel(-sigma)/dgumbel(u))))^2*k_1
                   k_7    <- -(y/mu)^(nu)*(tau-1)*(1-k_1)*((dgumbel(-sigma)*(1-exp(sigma))/dgumbel(u))-((tau-1)*dgumbel(-sigma)^2*(exp(-u)-1)/dgumbel(u)^2))
                   d2ldd2 <- k_3+k_4+k_5+k_6+k_7
                   return(d2ldd2)
                 },d2ldmdd = function(y,mu,sigma,nu,tau){ 
                   h_l     <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1     <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)*h_l
                   u       <- qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   l       <- (dgumbel(-sigma)*(tau-1))/dgumbel(u)
                   k_3     <- nu*(y^(nu)/mu^(nu+1))*(1+l)*(1-k_1)
                   k_4     <- nu*(y^(nu)/mu^(nu+1))*h_l*k_1*(1-(y/mu)^(nu)*(1+l))                   
                   d2ldmdd <- k_3-k_4
                   return(d2ldmdd)
                 },dldv = function(y,mu,sigma,nu,tau){  
                   h_l  <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1  <- exp(-(y/mu)^(nu)*h_l+sigma)
                   dldv <- (1/nu)+log(y)-log(mu)-(y/mu)^(nu)*log(y/mu)*h_l*(1-k_1)
                   return(dldv)
                 },dldt = function(y){
                   rep(0,length(y))
                 },d2ldv2 = function(y,mu,sigma,nu,tau){ 
                   h_l    <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1    <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2    <- (y/mu)^(nu)*log(y/mu)*h_l
                   d2ldv2 <- -(1/nu^2)-log(y/mu)*k_2*(1-k_1)-k_2^2*k_1
                   return(d2ldv2)
                 },d2ldt2 = function(y){
                   rep(0,length(y))
                 },d2ldmdv = function(y,mu,sigma,nu,tau){ 
                   h_l     <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1     <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)*h_l
                   t_1     <- -(1/mu)+k_2*(nu/mu)*(1-k_1)+(y/mu)^(nu)*(h_l/mu)*(1-k_1)
                   t_2     <- (nu/mu)*(y/mu)^(2*nu)*log(y/mu)*h_l^2*k_1 
                   d2ldmdv <- t_1+t_2
                   return(d2ldmdv)
                 },d2ldmdt = function(y){
                   rep(0,length(y))
                 },d2ldddv = function(y,mu,sigma,nu,tau){ 
                   h_l     <- sigma+qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   k_1     <- exp(-(y/mu)^(nu)*h_l+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)*h_l
                   u       <- qgumbel(tau*(1-pgumbel(-sigma))+pgumbel(-sigma))
                   l       <- (dgumbel(-sigma)*(tau-1))/dgumbel(u)
                   k_3     <- -(y/mu)^(nu)*log(y/mu)*(1+l)*(1-k_1)
                   k_4     <- k_1*k_2*(1-(y/mu)^(nu)*(1+l))
                   d2ldddv <- k_3+k_4
                   return(d2ldddv)
                 },d2ldddt = function(y){
                   rep(0, length(y))
                 },d2ldvdt = function(y){
                   rep(0, length(y))
                 },G.dev.incr = function(y,mu,sigma,nu,tau,...) -2*dRGGUMBEL(y,mu,sigma,nu,tau,log = TRUE),
                 rqres = expression(rqres(pfun = "pRGGUMBEL",type = "Continuous",y=y,mu=mu,sigma=sigma,nu=nu,tau=tau)), 
                 mu.initial = expression(mu <-rep(median(y),length(y))), sigma.initial = expression(sigma <- rep(0,length(y))), 
                 nu.initial = expression(nu<-rep(1,length(y))),tau.initial = expression(tau<-rep(0.5,length(y))), 
                 mu.valid = function(mu) all(mu>0), sigma.valid = function(sigma) TRUE, nu.valid = function(nu) all(nu>0),
                 tau.valid = function(tau) all((tau>0)&(tau<1)),y.valid = function(y) all(y>0)),class = c("gamlss.family","family"))}

