## This function called by robtest.R
jun.equation.2 <- function(beta, log.depth, y, phi, c.tukey) {
  return(equation.2(beta=beta, log.depth=log.depth, y=y, sig=phi, c.tukey=c.tukey))
}

## This function called by robtest.R
jun.all.vars <- function(y, log.depth, beta, phi, c.tukey.beta, c.tukey.phi) {
  return(all.vars(y=y, log.depth=log.depth, beta=beta, sig=phi, c.tukey.beta=c.tukey.beta, c.tukey.sig=c.tukey.phi))
}

## This function is called by robtest.R
jun.sig.rob.tukey <- function(phi, y, mu, c.tukey) {
  return(sig.rob.tukey(sig=phi, y=y, mu=mu, c.tukey=c.tukey))
}

## All functions following are for internal use only
varfunc <- function(mu,sig){mu+sig*mu^2}

loglkhd <- function(sig,y,mu){
  sum(lgamma(y+1/sig)-lgamma(1/sig)-lgamma(y+1)-(1/sig)*log(sig*mu+1)+y*log(sig*mu/(sig*mu+1)))
}

score.sig.ML <- function(sig,y,mu){
  sum(digamma(y+1/sig)-digamma(1/sig)-log(sig*mu+1)-sig*(y-mu)/(sig*mu+1))
}

info.sig.ML <- function(sig,y,mu){
  (-1/sig^2)*(sum(trigamma(y+1/sig))-length(y)*trigamma(1/sig))-sum((sig*mu^2+y)/(sig*mu+1)^2)
}

tukeypsi <- function(r,c.tukey){
  ifelse(abs(r)>c.tukey,0,((r/c.tukey)^2-1)^2*r)
}

E.tukeypsi.1 <- function(mui,sig,c.tukey){
  sqrtVmui <- sqrt(varfunc(mui,sig))
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
  }
}

E.tukeypsi.2 <- function(mui,sig,c.tukey){
  sqrtVmui <- sqrt(varfunc(mui,sig))
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)^2*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
  }
}

psi.sig.ML <- function(r,mu,sig){
  digamma(r*sqrt(mu*(sig*mu+1))+mu+1/sig)-sig*r*sqrt(mu/(sig*mu+1))-digamma(1/sig)-log(sig*mu+1)
}

ai.sig.tukey <- function(mui,sig,c.tukey){
  
  psi.sig.ML.mod <- function(j,mui,invsig){
    digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
  }
  sqrtVmui <- sqrt(mui*(sig*mui+1))
  invsig <- 1/sig
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))
  }
}

sig.rob.tukey <- function(sig,y,mu,c.tukey){
  r <- (y-mu)/sqrt(varfunc(mu,sig))
  wi <- tukeypsi(r=r,c.tukey=c.tukey)/r
  sum(wi*psi.sig.ML(r=r,mu=mu,sig=sig)-sapply(X=mu,FUN=ai.sig.tukey,sig=sig,c.tukey=c.tukey))
}

fullscore.sig <- function(y,mui,sigma){
  (digamma(y+1/sigma)-digamma(1/sigma)-log(sigma*mui+1)-sigma*(y-mui)/(sigma*mui+1))/(-sigma^2)
}

all.expectations.tukey <- function(mui,sigma,c.tukey.beta,c.tukey.sigma){
  expec <- list()
  sqrtVmui <- sqrt(varfunc(mui,sigma))
  j1.beta <- max(c(ceiling(mui-c.tukey.beta*sqrtVmui),0))
  j2.beta <- floor(mui+c.tukey.beta*sqrtVmui)
  if (j1.beta>j2.beta){
    expec$tukeypsi2 <- 0
    expec$psibetascoresig.beta <- 0
    expec$tukeypsi13 <- 0
    expec$psibetaminuspsisig <- 0
    j1.sigma <- max(c(ceiling(mui-c.tukey.sigma*sqrtVmui),0))
    j2.sigma <- floor(mui+c.tukey.sigma*sqrtVmui)
    if (j1.sigma>j2.sigma){
      expec$psibetascoresig.sigma <- 0
      expec$psiscoresig2 <- 0
      expec$psiscoresig13 <- 0
    } else {
      j12.sigma <- j1.sigma:j2.sigma
      probNB.sigma <- dnbinom(j12.sigma,mu=mui,size=1/sigma)
      resi.sigma <- (j12.sigma-mui)/sqrtVmui
      tukeyresi.sigma <- tukeypsi(r=resi.sigma,c.tukey=c.tukey.sigma)
      fullscoresig.sigma <- fullscore.sig(y=j12.sigma,mui=mui,sigma=sigma)
      expec$psibetascoresig.sigma <- sum(tukeyresi.sigma*fullscoresig.sigma*probNB.sigma)/sqrtVmui
      expec$psiscoresig2 <- sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma^2*probNB.sigma)
      expec$psiscoresig13 <- sum((tukeyresi.sigma/resi.sigma*fullscoresig.sigma)^2*probNB.sigma)+
        -sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)^2
    }
  } else {
    j12.beta <- j1.beta:j2.beta
    probNB.beta <- dnbinom(j12.beta,mu=mui,size=1/sigma)
    resi.beta <- (j12.beta-mui)/sqrtVmui
    tukeyresi.beta <- tukeypsi(r=resi.beta,c.tukey=c.tukey.beta)
    fullscoresig.beta <- fullscore.sig(y=j12.beta,mui=mui,sigma=sigma)
    expec$tukeypsi2 <- sum(tukeyresi.beta*(j12.beta-mui)*probNB.beta)/sqrtVmui^3
    expec$psibetascoresig.beta <- sum(tukeyresi.beta*fullscoresig.beta*probNB.beta)/sqrtVmui
    expec$tukeypsi13 <- (sum(tukeyresi.beta^2*probNB.beta)-sum(tukeyresi.beta*probNB.beta)^2)/sqrtVmui^2
    j1.sigma <- max(c(ceiling(mui-c.tukey.sigma*sqrtVmui),0))
    j2.sigma <- floor(mui+c.tukey.sigma*sqrtVmui)
    if (j1.sigma>j2.sigma){
      expec$psibetascoresig.sigma <- 0
      expec$psiscoresig2 <- 0
      expec$psibetaminuspsisig <- 0
      expec$psiscoresig13 <- 0
    } else {
      j12.sigma <- j1.sigma:j2.sigma
      probNB.sigma <- dnbinom(j12.sigma,mu=mui,size=1/sigma)
      resi.sigma <- (j12.sigma-mui)/sqrtVmui
      tukeyresi.sigma <- tukeypsi(r=resi.sigma,c.tukey=c.tukey.sigma)
      fullscoresig.sigma <- fullscore.sig(y=j12.sigma,mui=mui,sigma=sigma)
      expec$psibetascoresig.sigma <- sum(tukeyresi.sigma*fullscoresig.sigma*probNB.sigma)/sqrtVmui
      expec$psiscoresig2 <- sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma^2*probNB.sigma)
      if (j2.beta<j2.sigma){
        expec$psibetaminuspsisig <- (sum(tukeyresi.beta*tukeyresi.sigma[1:length(j12.beta)]/
                                           resi.sigma[1:length(j12.beta)]*fullscoresig.sigma[1:length(j12.beta)]*probNB.beta)+
                                       -sum(tukeyresi.beta*probNB.beta)*sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma))/sqrtVmui
      } else {
        expec$psibetaminuspsisig <- (sum(tukeyresi.beta[1:length(j12.sigma)]*tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)+
                                       -sum(tukeyresi.beta*probNB.beta)*sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma))/sqrtVmui
      }
      expec$psiscoresig13 <- sum((tukeyresi.sigma/resi.sigma*fullscoresig.sigma)^2*probNB.sigma)+
        -sum(tukeyresi.sigma/resi.sigma*fullscoresig.sigma*probNB.sigma)^2
    }
  }
  return(expec)
}

equation.2 <- function(beta, log.depth, y, sig, c.tukey){
  mu <- exp(beta + log.depth)
  r <- (y-mu) / sqrt(varfunc(mu,sig))
  res <- (tukeypsi(r=r,c.tukey=c.tukey) - sapply(X=mu, FUN=E.tukeypsi.1, sig=sig, c.tukey=c.tukey)) / sqrt(varfunc(mu,sig)) * mu * 1
  return(sum(res))
}

all.vars <- function(y, log.depth, beta, sig, c.tukey.beta, c.tukey.sig){
  
  n <- length(y)
  mu <- exp(beta + log.depth)
  expect <- sapply(X=mu, FUN=all.expectations.tukey, sigma=sig, c.tukey.beta=c.tukey.beta, c.tukey.sigma=c.tukey.sig)
  M11 <- sum(as.numeric(unlist(expect['tukeypsi2',])*mu^2)) / n
  M12 <- sum(as.numeric(unlist(expect['psibetascoresig.beta',])*mu)) / n
  M21 <- sum(as.numeric(unlist(expect['psibetascoresig.sigma',])*mu)) / n
  M22 <- sum(as.numeric(unlist(expect['psiscoresig2',]))) / n
  
  Q11 <- sum(as.numeric(unlist(expect['tukeypsi13',])*mu^2)) / n
  Q12 <- sum(as.numeric(unlist(expect['psibetaminuspsisig',])*mu)) / n
  Q22 <- sum(as.numeric(unlist(expect['psiscoresig13',]))) / n
  
  fullM <- rbind(cbind(M11,M12),t(c(M21,M22)))
  fullQ <- rbind(cbind(Q11,Q12),t(c(Q12,Q22)))
  vars <- solve(fullM)%*%fullQ%*%solve(fullM) / n
  
  return(vars)
}
