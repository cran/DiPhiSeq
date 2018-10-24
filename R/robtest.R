#' Calls the robnb function to estimate the coefficients,
#' and then construct the statistical tests for DD and DE.
#' It works for a single gene. y1 and y2 are count vectors for a single gene.
#' diphiseq calls this function to do the calculation for each gene.
#' Normal users often don't need to use this function directly.
#'
#' @param y1 counts from group 1. a vector.
#' @param log.depth1 log(sequencing depths) for samples in group 1. a vector.
#' @param y2 counts from group 2. a vector.
#' @param log.depth2 log(sequencing depths) for samples in group 2. a vector.
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
#' res <- robtest(y1, log(d1), y2, log(d2))
#' @export
robtest <- function(y1, log.depth1, y2, log.depth2, c.tukey.beta=4, c.tukey.phi=4) {
  res1 <- robnb(y=y1, log.depth=log.depth1, c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi)
  res2 <- robnb(y=y2, log.depth=log.depth2, c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi)

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
#' @param log.depth Vector of log(sequencing depths).
#' @param c.tukey.beta The c value for beta in Huber function. The default value should be appropriate
#'   for most datasets.
#' @param c.tukey.phi The c value for phi in Huber function. The default value should be appropriate
#'   for most datasets.
#' @param phi.ini The initial value of phi.
#' @param alpha A positive value for setting initial values. The default value is usually appropriate.
#' @param minphi A searching parameter for Algorithm 1 (check the algorithm for details.)
#'   The default value is usually appropriate.
#' @param maxphi A searching parameter for Algorithm 1 (check the algorithm for details.)
#'   The default value is usually appropriate.
#' @param maxit Maximum number of iterations for the outer loop.
#'   The default value is usually appropriate.
#' @param maxit.beta Maximum number of iterations for the inner loop of solving beta.
#'   The default value is usually appropriate.
#' @param maxit.phi Maximum number of iterations for the inner loop of solving phi.
#'   The default value is usually appropriate.
#' @param tol.beta The numerical tolerance of solving beta.
#'   The default value is usually appropriate.
#' @param tol.phi The numerical tolerance of solving phi.
#'   The default value is usually appropriate.
#' @return A list that contains the elements:
#'   \code{beta}: the estimated (log) expression.
#'   \code{phi}: the estimated dispersion.
#'   \code{fconv}: flag of the convergence of the search.
#'   \code{vars}: the variance-covariance matrix of the estimates.
#'   \code{sd.beta}: the standard error of beta.
#'   \code{sd.phi}: the standard error of phi.
#'   \code{y}: the input y value.
#'   \code{log.depth}: log(sequencing depth).
#' @examples
#' d <- runif(10, min=1, max=2)
#' y <- rnbinom(10, size=1, mu=d*50)
#' res <- robnb(y, log(d))
#' @export
robnb <- function(y, log.depth, c.tukey.beta=4, c.tukey.phi=4, phi.ini=0.5,
                            alpha=0.2, minphi=0.01, maxphi=5, maxit=30, maxit.beta=30, maxit.phi=30,
                            tol.beta=0.01, tol.phi=0.005)
{
  beta.ini <- ini.value.beta(y, log.depth, alpha=alpha)

  rob.sol <- rob.glm(y=y, log.depth=log.depth, beta=beta.ini, phi=phi.ini, c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi,
                         minphi=minphi, maxphi=maxphi, maxit=maxit, maxit.beta=maxit.beta, maxit.phi=maxit.phi,
                         tol.beta=tol.beta, tol.phi=tol.phi)
  beta <- rob.sol$beta
  phi <- rob.sol$phi
  fconv <- rob.sol$fconv

  vars <- jun.all.vars(y=y, log.depth=log.depth, beta=beta, phi=phi, c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi)
  sd.beta <- sqrt(vars[1, 1])
  sd.phi <- sqrt(vars[2, 2])

  return(list(beta=beta, phi=phi, fconv=fconv, vars=vars, sd.beta=sd.beta, sd.phi=sd.phi, y=y, log.depth=log.depth))
}

## search an interval of phi that contains the next solution
set.minphi.maxphi <- function(phi, log.depth, y, beta, c.tukey.phi, gint=1.2, nint=40) {
  sign0 <- sign(jun.sig.rob.tukey(phi=phi, y=y, mu=exp(beta + log.depth), c.tukey=c.tukey.phi))
  for (i in 1 : nint) {
    sign1 <- sign(jun.sig.rob.tukey(phi=phi / gint ^ i, y=y, mu=exp(beta + log.depth), c.tukey=c.tukey.phi))
    if (sign1 * sign0 <= 0) {minphi <- phi / gint ^ i; maxphi <- phi; break}
    
    sign2 <- sign(jun.sig.rob.tukey(phi=phi * gint ^ i, y=y, mu=exp(beta + log.depth), c.tukey=c.tukey.phi))
    if (sign2 * sign0 <= 0) {minphi <- phi; maxphi <- phi * gint ^ i; break}
  }
  if (i == nint) {stop("Cannot identify the correct range of phi!")}
  
  return(c(minphi, maxphi))
}

## search an interval of beta that contains the next solution
set.minbeta.maxbeta <- function(beta, log.depth, y, phi, c.tukey.beta, gint=0.1, nint=40) {
  sign0 <- sign(jun.equation.2(beta, log.depth, y, phi, c.tukey.beta))
  for (i in 1 : nint) {
    sign1 <- sign(jun.equation.2(beta - gint * i, log.depth, y, phi, c.tukey.beta))
    if (sign1 * sign0 <= 0) {minbeta <- beta - gint * i; maxbeta <- beta; break}
    
    sign2 <- sign(jun.equation.2(beta + gint * i, log.depth, y, phi, c.tukey.beta))
    if (sign2 * sign0 <= 0) {minbeta <- beta; maxbeta <- beta + gint * i; break}
  }
  if (i == nint) {stop("Cannot identify the correct range of beta!")}
  
  return(c(minbeta, maxbeta))
}

## given initial values of beta and phi, solve the glm problem
rob.glm <- function(y, log.depth, beta, phi, c.tukey.beta, c.tukey.phi,
                    minphi, maxphi, maxit, maxit.beta, maxit.phi, tol.beta, tol.phi) {
  fconv <- FALSE
  for (it in 1 : maxit) {
    cat("\n Iteration ", it, " ")
    
    phi0 <- phi
    beta0 <- beta
    
    phi.interval <- set.minphi.maxphi(phi=phi, log.depth=log.depth, y=y, beta=beta, c.tukey.phi=c.tukey.phi, gint=1.2, nint=40)
    phi <- uniroot(f=jun.sig.rob.tukey, interval=phi.interval, tol=tol.phi, maxiter=maxit.phi, mu=exp(beta + log.depth), y=y, c.tukey=c.tukey.phi)$root
    if (is.na(phi) | phi>maxphi) {phi <- maxphi}
    cat(".")
    
    beta.interval <- set.minbeta.maxbeta(beta=beta, log.depth=log.depth, y=y, phi=phi, c.tukey.beta=c.tukey.beta, gint=0.1, nint=40)
    beta <- uniroot(f=jun.equation.2, interval=beta.interval, tol=tol.beta, maxiter=maxit.beta, log.depth=log.depth, y=y, phi=phi, c.tukey=c.tukey.beta)$root
    cat(".", fill=TRUE)
    
    cat("phi =", phi, fill=TRUE)
    cat("beta =", beta, fill=TRUE)
    
    if (abs(phi - phi0) < tol.phi & abs(beta - beta0) < tol.beta) {fconv <- TRUE; break}
  }
  
  return(list(beta=beta, phi=phi, fconv=fconv))
}

## this function gives a high-breakpoint naive estimate of the initial value of beta
ini.value.beta <- function(y, log.depth, alpha=0.2) {
  ord <- order(y / exp(log.depth))
  n <- length(y)
  ndel <- floor(n * alpha)
  if (ndel >= floor(n / 2)) {ndel <- floor(n / 2) - 1}
  if (ndel < 0) {ndel <- 0}
  tokeep <- ord[(ndel + 1) : (n - ndel)]
  beta.ini <- log(max(sum(y[tokeep]) / sum(exp(log.depth)[tokeep]), 0.01))
  return(beta.ini)
}

## this function gives a high-breakpoint naive estimate of the initial value of sigma
## We do not directly use this estimate as it may be too unstable. Instead,
## We use the median value of all genes
ini.value.phi.each.gene <- function(y, log.depth, beta.ini, alpha=0.2) {
  mu <- exp(log.depth + beta.ini)
  phi.s <- ((y - mu) ^ 2 - mu) / (mu ^ 2)
  phi.each.gene <- mean(phi.s, trim=alpha, na.rm=TRUE)
  return(phi.each.gene)
}

## this function gives the final initial value of phi
## here two new parameters are introduced: min and max start value for the phi
ini.value.phi <- function(countmat, classlab, log.depth, alpha=0.2) {
  phi.each.gene <- matrix(NA, nrow=nrow(countmat), ncol=2)
  for (i in 1 : nrow(countmat)) {
    for (j in 1 : 2) {
      beta.ini <- ini.value.beta(y=countmat[i, classlab == j], log.depth=log.depth[classlab == j], alpha=alpha)
      phi.each.gene[i, j] <- ini.value.phi.each.gene(y=countmat[i, classlab == j], 
          log.depth=log.depth[classlab == j], beta.ini=beta.ini, alpha=alpha)
    }
  }
  phi.ini <- median(phi.each.gene)

  return(phi.ini)
}
