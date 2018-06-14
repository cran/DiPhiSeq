#' Calls the robnb function to estimate the coefficients,
#' and then construct the statistical tests for DD and DE.
#' It works for a single gene. y1 and y2 are count vectors for a single gene.
#' diphiseq calls this function to do the calculation for each gene.
#' Normal users often don't need to use this function directly.
#'
#' @param y1 counts from group 1. a vector.
#' @param d1 sequencing depth for samples in group 1. a vector.
#' @param y2 counts from group 2. a vector.
#' @param d2 sequencing depth for samples in group 2. a vector.
#' @param c.tukey.beta The c value for beta in Huber function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @param c.tukey.phi The c value for phi in Huber function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @return A vector that contains the elements:
#'   \code{phi1}: the estimated dispersion of sample group 1.
#'   \code{phi2}: the estimated dispersion of sample group 2.
#'   \code{beta1}: the estimated (log) expression of sample group 1.
#'   \code{beta2}: the estimated (log) expression of sample group 2.
#'   \code{statistic.phi}: the z statistic for DD.
#'   \code{statistic.beta}: the z statistic for DE.
#'   \code{p.value.phi}: the p value for DD.
#'   \code{p.value.beta}: the p value for DE.
#' @examples
#' d1 <- runif(10, min=1, max=2)
#' d2 <- runif(15, min=1, max=2)
#' y1 <- rnbinom(10, size=1, mu=d1*50)
#' y2 <- rnbinom(15, size=1, mu=d2*50)
#' res <- robtest(y1, d1, y2, d2)
robtest <- function(y1, d1, y2, d2, c.tukey.beta=4, c.tukey.phi=4)
{
  c.tukey.sig <- c.tukey.phi # in robnb, c.,tukey.sig is used.

  res1 <- robnb(y1, d1, c.tukey.beta=c.tukey.beta, c.tukey.sig=c.tukey.sig)
  res2 <- robnb(y2, d2, c.tukey.beta=c.tukey.beta, c.tukey.sig=c.tukey.sig)

  statistic.phi <- (res2$phi - res1$phi) / sqrt(res2$sd.phi^2 + res1$sd.phi^2)
  statistic.beta <- (res2$beta - res1$beta) / sqrt(res2$sd.beta^2 + res1$sd.beta^2)

  p.value.phi <- 2 * pnorm(-abs(statistic.phi))
  p.value.beta <- 2 * pnorm(-abs(statistic.beta))

  res <- c(phi1=res1$phi, phi2=res2$phi, beta1=res1$beta, beta2=res2$beta,
           statistic.phi=statistic.phi, statistic.beta=statistic.beta,
           p.value.phi=p.value.phi, p.value.beta=p.value.beta)
}

#' Calculates the estimate and standard error of beta and phi.
#' It takes as input counts from one group of samples for a single gene.
#' This function is the core underlining function of the whole package.
#' A significant part of the code is edited based on William H. Aeberhard's glmrob.nb R function;
#' we appreciate them very much for sharing their code online.
#' This function also implement Algorithm 1 of our submitted paper about DiPhiSeq.
#' This function is called by robtest.
#' Most users don't need to call this function directly.
#'
#' @param y A count vector.
#' @param d Vector of sequencing depths.
#' @param c.tukey.beta The c value for beta in Huber function. The default value should be appropriate
#'   for most datasets.
#' @param c.tukey.sig The c value for phi in Huber function. The default value should be appropriate
#'   for most datasets.
#' @param alpha A positive value for setting initial values. The default value is usually appropriate.
#' @param minsig A searching parameter for Algorithm 1 (check the algorithm for details.)
#'   The default value is usually appropriate.
#' @param maxsig A searching parameter for Algorithm 1 (check the algorithm for details.)
#'   The default value is usually appropriate.
#' @param maxit Maximum number of iterations for the outer loop.
#'   The default value is usually appropriate.
#' @param maxit.beta Maximum number of iterations for the inner loop of solving beta.
#'   The default value is usually appropriate.
#' @param maxit.sig Maximum number of iterations for the inner loop of solving phi (named sig in this function).
#'   The default value is usually appropriate.
#' @param tol.beta The numerical tolerance of solving beta.
#'   The default value is usually appropriate.
#' @param tol.sig The numerical tolerance of solving phi (named sig in this function).
#'   The default value is usually appropriate.
#' @param sig.ini A searching parameter for Algorithm 1 (check the algorithm for details.)
#'   The default value is usually appropriate.
#' @return A list that contains the elements:
#'   \code{beta}: the estimated (log) expression.
#'   \code{phi}: the estimated dispersion.
#'   \code{fconv}: flag of the convergence of the search.
#'   \code{vars}: the variance-covariance matrix of the estimates.
#'   \code{sd.beta}: the standard error of beta.
#'   \code{sd.phi}: the standard error of phi.
#'   \code{y}: the input y value.
#'   \code{d}: the input d value.
#'   \code{D}: log(d).
#' @examples
#' d <- runif(10, min=1, max=2)
#' y <- rnbinom(10, size=1, mu=d*50)
#' res <- robnb(y, d)
robnb <- function(y, d, c.tukey.beta=4, c.tukey.sig=4,
                            alpha=0.2, minsig=0.01, maxsig=5, maxit=30, maxit.beta=30, maxit.sig=30,
                            tol.beta=0.01, tol.sig=0.005, sig.ini=0.5)
{

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

	equation.2 <- function(beta, D, y, sig, c.tukey){
		mu <- exp(beta + D)
		r <- (y-mu) / sqrt(varfunc(mu,sig))
		res <- (tukeypsi(r=r,c.tukey=c.tukey) - sapply(X=mu, FUN=E.tukeypsi.1, sig=sig, c.tukey=c.tukey)) / sqrt(varfunc(mu,sig)) * mu * 1
		return(sum(res))
	}

	set.minsig.maxsig <- function(sig, D, y, beta, c.tukey.sig, gint=1.2, nint=40){
		sign0 <- sign(sig.rob.tukey(sig=sig, y=y, mu=exp(beta + D), c.tukey=c.tukey.sig))
		for (i in 1 : nint)
		{
			sign1 <- sign(sig.rob.tukey(sig=sig / gint ^ i, y=y, mu=exp(beta + D), c.tukey=c.tukey.sig))
			if (sign1 * sign0 <= 0) {minsig <- sig / gint ^ i; maxsig <- sig; break}

			sign2 <- sign(sig.rob.tukey(sig=sig * gint ^ i, y=y, mu=exp(beta + D), c.tukey=c.tukey.sig))
			if (sign2 * sign0 <= 0) {minsig <- sig; maxsig <- sig * gint ^ i; break}
		}
		if (i == nint) {stop("Cannot identify the correct range of sigma!")}

		return(c(minsig, maxsig))
	}

	set.minbeta.maxbeta <- function(beta, D, y, sig, c.tukey.beta, gint=0.1, nint=40){
		sign0 <- sign(equation.2(beta, D, y, sig, c.tukey.beta))
		for (i in 1 : nint)
		{
			sign1 <- sign(equation.2(beta - gint * i, D, y, sig, c.tukey.beta))
			if (sign1 * sign0 <= 0) {minbeta <- beta - gint * i; maxbeta <- beta; break}

			sign2 <- sign(equation.2(beta + gint * i, D, y, sig, c.tukey.beta))
			if (sign2 * sign0 <= 0) {minbeta <- beta; maxbeta <- beta + gint * i; break}
		}
		if (i == nint) {stop("Cannot identify the correct range of beta!")}

		return(c(minbeta, maxbeta))
	}

	rob.glm <- function(y, D, beta, sig, c.tukey.beta, c.tukey.sig,
                              minsig, maxsig, maxit, maxit.beta, maxit.sig, tol.beta, tol.sig){
		fconv <- FALSE
		for (it in 1 : maxit)
		{
			cat("\n Iteration ", it, " ")

			sig0 <- sig
			beta0 <- beta

			sig.interval <- set.minsig.maxsig(sig=sig, D=D, y=y, beta=beta, c.tukey.sig=c.tukey.sig, gint=1.2, nint=40)
			sig <- uniroot(f=sig.rob.tukey, interval=sig.interval, tol=tol.sig, maxiter=maxit.sig, mu=exp(beta + D), y=y, c.tukey=c.tukey.sig)$root
			if (is.na(sig) | sig>maxsig) {sig <- maxsig}
			cat(".")

			beta.interval <- set.minbeta.maxbeta(beta=beta, D=D, y=y, sig=sig, c.tukey.beta=c.tukey.beta, gint=0.1, nint=40)
			beta <- uniroot(f=equation.2, interval=beta.interval, tol=tol.beta, maxiter=maxit.beta, D=D, y=y, sig=sig, c.tukey=c.tukey.beta)$root
			cat(".", fill=TRUE)

			cat("sig =", sig, fill=TRUE)
			cat("beta =", beta, fill=TRUE)

			if (abs(sig - sig0) < tol.sig & abs(beta - beta0) < tol.beta) {fconv <- TRUE; break}
		}

		return(list(beta=beta, sig=sig, fconv=fconv))
	}

	all.vars <- function(y, D, beta, sig, c.tukey.beta, c.tukey.sig){

		n <- length(y)
		mu <- exp(beta + D)
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

	ini.values <- function(y, D, alpha=0.2){
		ord <- order(y / exp(D))
		n <- length(y)
		ndel <- floor(n * alpha)
		if (ndel >= floor(n / 2)) {ndel <- floor(n / 2) - 1}
		if (ndel < 0) {ndel <- 0}
		tokeep <- ord[(ndel + 1) : (n - ndel)]
		beta.ini <- log(max(sum(y[tokeep]) / sum(exp(D)[tokeep]), 0.01))
		return(beta.ini)
	}

  #d <- exp(log(d) - mean(log(d))) # this is a line I commented from Alicia's code
  D <- log(d)

  beta.ini <- ini.values(y, D, alpha=alpha)

  rob.sol <- rob.glm(y=y, D=D, beta=beta.ini, sig=sig.ini, c.tukey.beta=c.tukey.beta, c.tukey.sig=c.tukey.sig,
                         minsig=minsig, maxsig=maxsig, maxit=maxit, maxit.beta=maxit.beta, maxit.sig=maxit.sig,
                         tol.beta=tol.beta, tol.sig=tol.sig)
  beta <- rob.sol$beta
  sig <- rob.sol$sig
  fconv <- rob.sol$fconv

  vars <- all.vars(y=y, D=D, beta=beta, sig=sig, c.tukey.beta=c.tukey.beta, c.tukey.sig=c.tukey.sig)
  sd.beta <- sqrt(vars[1, 1])
  sd.sig <- sqrt(vars[2, 2])

  return(list(beta=beta, phi=sig, fconv=fconv, vars=vars, sd.beta=sd.beta, sd.phi=sd.sig, y=y, d=d, D=D))
}
