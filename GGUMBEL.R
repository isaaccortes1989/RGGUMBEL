#===========================================#
#===Install and load the following packages
#===========================================#
library(evd)
library(gamlss)

dGGUMBEL = function(x,mu=1,sigma=1,nu=1,log=FALSE){
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(mu<=0)) stop(paste("mu must be positive","\n",""))
lfy <- log(nu)+(nu-1)*log(x)-nu*log(mu)-log(1-pgumbel(-sigma))+dgumbel((x/mu)^(nu)-sigma,log = TRUE)
  if(log== FALSE)
    fy <- exp(lfy)
  else fy <- lfy
  return(fy)
}

pGGUMBEL = function(q,mu=1,sigma=1,nu=1,lower.tail=TRUE,log.p=FALSE){
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(mu<=0))  stop(paste("mu must be positive","\n",""))
lGy     <- log(pgumbel((q/mu)^(nu)-sigma)-pgumbel(-sigma))-log(1-pgumbel(-sigma)) 
  if(log.p == FALSE)
    cdf <- exp(lGy)  
  else cdf <- lGy
  if(lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1-cdf
  return(cdf)
}

qGGUMBEL = function(p,mu=1,sigma=1,nu=1,lower.tail=TRUE,log.p=FALSE){
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(mu<=0))  stop(paste("mu must be positive","\n",""))
  if(log.p == TRUE)
    p <- exp(p)  
  else p <- p
  if(lower.tail == TRUE)
    p <- p
  else p <- 1-p
  if(any(p<0)|any(p>1)) stop(paste("p must be between 0 and 1","\n",""))
  Q       <- mu*(sigma+qgumbel(p*(1-pgumbel(-sigma))+pgumbel(-sigma)))^(1/nu)
  return(Q)
}

rGGUMBEL = function(n,mu=1,sigma=1,nu=1){
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(mu<=0))  stop(paste("mu must be positive","\n",""))
  U <- runif(n)
  X <- mu*(sigma+qgumbel(U*(1-pgumbel(-sigma))+pgumbel(-sigma)))^(1/nu)
  return(X)
} 


GGUMBEL<-function(mu.link = "identity", sigma.link = "identity",nu.link = "identity"){
  mstats <- checklink("mu.link","GGUMBEL",substitute(mu.link),c("inverse","log","identity","own"))
  dstats <- checklink("sigma.link","GGUMBEL",substitute(sigma.link),c("inverse","log","identity","own"))
  vstats <- checklink("nu.link","GGUMBEL",substitute(nu.link),c("log","identity"))
  structure(list(family = c("GGUMBEL","Generalized Gumbel"), 
                 parameters = list(mu=TRUE,sigma=TRUE,nu=TRUE),nopar=3,type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)), 
                 mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun, 
                 mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv, 
                 mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta,
                 dldm = function(y,mu,sigma,nu){
                   k_1  <- exp(-(y/mu)^(nu)+sigma)
                   dldm <- -(nu/mu)+nu*(y^(nu)/mu^(nu+1))*(1-k_1)
                   return(dldm)
                 },d2ldm2 = function(y,mu,sigma,nu){
                   k_1    <- exp(-(y/mu)^(nu)+sigma)
                   d2ldm2 <-(nu/mu^2)-nu*(nu+1)*(y^(nu)/mu^(nu+2))*(1-k_1)-nu^2*(y^(2*nu)/mu^(2*nu+2))*k_1
                   return(d2ldm2)
                 },dldd = function(y,mu,sigma,nu){
                   k_1  <- exp(-(y/mu)^(nu)+sigma)
                   dldd <- (dgumbel(-sigma)/(1-pgumbel(-sigma)))+1-k_1
                   return(dldd)
                 },d2ldd2 = function(y,mu,sigma,nu){
                   k_1  <- exp(-(y/mu)^(nu)+sigma)
                   k_2  <- dgumbel(-sigma)/(1-pgumbel(-sigma))
                   d2ldd2 <- -k_2^2+k_2*(exp(-sigma)-1)-k_1
                   return(d2ldd2)
                 },d2ldmdd = function(y,mu,sigma,nu){
                   k_1     <- exp(-(y/mu)^(nu)+sigma)
                   d2ldmdd <- -nu*(y^(nu)/mu^(nu+1))*k_1
                   return(d2ldmdd)
                 },dldv = function(y,mu,sigma,nu){
                   k_1  <- exp(-(y/mu)^(nu)+sigma)
                   dldv <- (1/nu)+log(y)-log(mu)+(y/mu)^(nu)*log(y/mu)*(k_1-1)
                   return(dldv)
                 },d2ldv2 = function(y,mu,sigma,nu){
                   k_1    <- exp(-(y/mu)^(nu)+sigma)
                   k_2    <- (y/mu)^(nu)*log(y/mu)
                   d2ldv2 <- -(1/nu^2)+log(y/mu)*k_2*(k_1-1)-k_2^(2)*k_1
                   return(d2ldv2)
                 },d2ldmdv = function(y,mu,sigma,nu){
                   k_1     <- exp(-(y/mu)^(nu)+sigma)
                   k_2     <- y^(nu)/mu^(nu+1)
                   t_1     <- -(1/mu)+k_2*(1-k_1)+nu*log(y)*k_2*(1-k_1)
                   t_2     <- -nu*log(mu)*k_2*(1-k_1)+nu*k_2*k_1*(y/mu)^(nu)*log(y/mu) 
                   d2ldmdv <- t_1+t_2
                   return(d2ldmdv)
                 },d2ldddv = function(y,mu,sigma,nu){
                   k_1     <- exp(-(y/mu)^(nu)+sigma)
                   k_2     <- (y/mu)^(nu)*log(y/mu)
                   d2ldddv <- k_1*k_2
                   return(d2ldddv)
                 },G.dev.incr = function(y,mu,sigma,nu,...) -2*dGGUMBEL(y,mu,sigma,nu,log = TRUE),
                 rqres = expression(rqres(pfun = "pGGUMBEL",type = "Continuous",y=y,mu=mu,sigma=sigma,nu=nu)), 
                 mu.initial = expression(mu <-rep(sqrt(mean(y^2)),length(y))), sigma.initial = expression(sigma <- rep(0,length(y))), 
                 nu.initial = expression(nu<-rep(1,length(y))), 
                 mu.valid = function(mu) all(mu>0), sigma.valid = function(sigma) TRUE, nu.valid = function(nu) all(nu>0),
                 y.valid = function(y) all(y>0)),class = c("gamlss.family","family"))}



