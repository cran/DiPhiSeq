#' Main function. For most users, this function is all what they
#' need for the analysis.
#'
#' @param countmat A count matrix. Rows are genes, and columns are samples. Each element/count is
#'   the number of reads mapped to a gene in a sample.
#' @param classlab The class labels. A vector whose elements are of value 1 or 2.
#' @param depth Sequencing depth. A vector of sequencing depth. Users are encouraged to
#'   provide estimated values from edgeR, DESeq, PoissonSeq, or other software that they prefer.
#'   If no values are provided, depth estimated by total counts will be used.
#' @param c.tukey.beta The c value for beta in Tukey's biweight function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @param c.tukey.phi The c value for phi in Tukey's biweight function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @param phi.ini The initial value for phi to start search. If phi.in == 'adaptive' (the default value), 
#'   the algorithm will adaptively
#'   choose a value (check Algorithm 1 for details). Otherwise, a positive numeric value (such as 0.5)
#'   should be given.
#' @return A List that contains the following elements:
#' \code{tab}: This is a data frame that contains the main results of this package. It has the following columns:
#'   \code{phi1}: the estimated dispersion of sample group 1.
#'   \code{phi2}: the estimated dispersion of sample group 2.
#'   \code{beta1}: the estimated (log) expression of sample group 1.
#'   \code{beta2}: the estimated (log) expression of sample group 2.
#'   \code{statistic.phi}: the z statistic for DD.
#'   \code{statistic.beta}: the z statistic for DE.
#'   \code{p.value.phi}: the p value for DD.
#'   \code{p.value.beta}: the p value for DE.
#'   \code{fdr.phi}: the FDR for DD.
#'   \code{fdr.beta}: the FDR for DE.
#' \code{log.depth}: A vector of the log sequencing depths.
#' \code{countmat}: This is the count matrix that the user provided.
#' \code{classlab}: This is the vecgtor of class labels that the user provided.
#' \code{phi.ini}: The initial searching value of the dispersion parameter.
#' \code{mumat}: This is the (estimated) mu (expected expression) matrix, of the same size as countmat.
#'   In another word, E(countmat)=mumat.
#' \code{phimat}: This is the (estimated) phi matrix, of the same size as countmat.
#'   In another word, counmat ~ negative binomial(mumat, phimat).
#' @examples
#' countmat <- matrix(rnbinom(100, size=1, mu=50), nrow=4, ncol=25)
#' classlab <- c(rep(1, 10), rep(2, 15))
#' res <- diphiseq(countmat, classlab)
#' 
#' countmat <- matrix(rnbinom(100, size=1, mu=50), nrow=4, ncol=25)
#' classlab <- c(rep(1, 10), rep(2, 15))
#' res <- diphiseq(countmat, classlab, phi.ini=0.5)
#' @export
diphiseq <- function(countmat, classlab, depth=NULL, c.tukey.beta=4, c.tukey.phi=4, phi.ini='adaptive')
{
  ## check arguments
  if (ncol(countmat) != length(classlab)){
    stop("The dimensions of countmat and classlab do not agree!")
  }

  if (sum(countmat < 0) > 0){
    stop("There are negative values in the countmat: please note that diphiseq only accepts original counts,
         not normalized or log transformed data!")
  }

  if (sum(classlab == 1) + sum(classlab == 2) < length(classlab)){
    stop("Values in classlab can only be 1 or 2!")
  }

  if (sum(classlab == 1) == 0 | sum(classlab == 2) == 0){
    stop("classlab must contain both 1 and 2!")
  }

  if (!is.null(depth)){
    if (length(depth) != length(classlab)){
      stop("lengths of depth and classlab do not agree!")
    }
    if (min(depth) <= 0){
      stop("depth must all be positive values")
    }
  }

  if (c.tukey.beta <= 0){
    stop("c.tukey.beta must be positive!")
  }
  if (c.tukey.phi <= 0){
    stop("c.tukey.phi must be positive!")
  }
  
  ## sequencing depth
  if (is.null(depth)){
    cat("Sequencing depths not given. Estimated by total counts.\n")
    depth <- colMeans(countmat)
  } else if (sum(depth <= 0) > 0) {
    stop("Sequencing depths must be positive values!")
  }
  depth <- exp(log(depth) - mean(log(depth))) # adjust sequencing depth so that they have geometric mean 0
  log.depth <- log(depth)

  ## initial value of phi
  if (phi.ini != 'adaptive') {
    if (!is.numeric(phi.ini)) {
      stop("phi.ini must be 'adaptive' or a positive numeric value!")
    } else if (phi.ini < 0) {
      stop("phi.ini must be 'adaptive' or a positive numeric value!")
    }
  }
  
  ## phi.ini
  if (phi.ini == 'adaptive') {
    phi.ini <- ini.value.phi(countmat=countmat, classlab=classlab, log.depth=log.depth)
    cat("The adaptive initial value of phi for this dataset is", phi.ini, fill=TRUE)
  }
  
  ## adjust back if it is out of range 0.1 and 1.0
  min.st.phi <- 0.1
  max.st.phi <- 1
  if (phi.ini < min.st.phi | phi.ini > max.st.phi) {
    if (phi.ini < min.st.phi) {phi.ini <- min.st.phi}
    if (phi.ini > max.st.phi) {phi.ini <- max.st.phi}
    cat("The adaptive initial value of phi is suggested to be between", min.st.phi, "and", max.st.phi, "; so", phi.ini, 
        "will be used. Typically, the results of DiPhiSeq is insensitive to the starting values anyway.", fill=TRUE)
  }

  ## start computing
  res <- matrix(NA, nrow=nrow(countmat), ncol=8)
  for (i in 1 : nrow(countmat))
  {
    tryCatch({
      res[i, ] <- robtest(y1=countmat[i, classlab == 1], log.depth1=log.depth[classlab == 1],
                          y2=countmat[i, classlab == 2], log.depth2=log.depth[classlab == 2],
                          c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi)
      if (i %% 10 == 0) {
        cat(i, "genes processed.\n")
      }
    }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  }

  ## convert the results into a data frame
  res <- data.frame(res)
  colnames(res) <- c('phi1', 'phi2', 'beta1', 'beta2', 'statistic.phi', 'statistic.beta',
                    'p.value.phi', 'p.value.beta')

  ## convert to fdr values
  res$fdr.phi <- p.adjust(res$p.value.phi, method="fdr")

  # FDR corrected p-values for DE test
  res$fdr.beta <- p.adjust(res$p.value.beta, method="fdr")
  
  ## calculate mumat and phimat
  mumat <- phimat <- countmat
  for (i in 1 : nrow(mumat)) {
    mumat[i, classlab==1] <- exp(res$beta1[i] + log.depth[classlab==1])
    mumat[i, classlab==2] <- exp(res$beta2[i] + log.depth[classlab==2])
    phimat[i, classlab==1] <- res$phi1[i]
    phimat[i, classlab==2] <- res$phi2[i]
  }
  
  res.list <- list(tab=res, log.depth=log.depth, countmat=countmat, classlab=classlab, phi.ini=phi.ini, 
                   mumat=mumat, phimat=phimat)
  return(res.list)
}


#' Give a rough estimate of the proportion of outliers in the data based on the results of DiPhiSeq.
#'
#' @param diphiseq.res The results given by running diphiseq.
#' @param fdr.cutoff The cutoff for FDR. 
#' @return a numeric value. The estimated proportion of outliers under the FDR cutoff in the data.
#' @examples
#' countmat <- matrix(rnbinom(100, size=1, mu=50), nrow=4, ncol=25)
#' classlab <- c(rep(1, 10), rep(2, 15))
#' res <- diphiseq(countmat, classlab)
#' outlier.proportion <- outprop(res)
#' @export
outprop <- function(diphiseq.res, fdr.cutoff=0.1) {
  varmat <- diphiseq.res$mumat + diphiseq.res$mumat ^ 2 * diphiseq.res$phimat
  rmat <- (diphiseq.res$countmat - diphiseq.res$mumat) / sqrt(varmat)
  rvec <- c(rmat)
  pvals <- pnorm(-abs(rvec)) * 2
  fdrs <- p.adjust(pvals, method='fdr')
  outlier.prop <- mean(fdrs < fdr.cutoff)
  return(outlier.prop)
}
