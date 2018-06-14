#' Main function. For most users, this function is all what they
#' need for the analysis.
#'
#' @param countmat A count matrix. Rows are genes, and columns are samples. Each element/count is
#'   the number of reads mapped to a gene in a sample.
#' @param classlab The class labels. A vector whose elements are of value 1 or 2.
#' @param depth Sequencing depth. A vector of sequencing depth. Users are encouraged to
#'   provide estimated values from edgeR, DESeq, PoissonSeq, or other software that they prefer.
#'   If no values are provided, depth estimated by total counts will be used.
#' @param c.tukey.beta The c value for beta in Huber function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @param c.tukey.phi The c value for phi in Huber function. The default value, 4, is typically
#'   regarded as appropriate and should work for most datasets.
#' @return A data.frame that contains the following columns:
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
#' @examples
#' countmat <- matrix(rnbinom(100, size=1, mu=50), nrow=4, ncol=25)
#' classlab <- c(rep(1, 10), rep(2, 15))
#' res <- diphiseq(countmat, classlab)
diphiseq <- function(countmat, classlab, depth=NULL, c.tukey.beta=4, c.tukey.phi=4)
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
    stop("c.tukey.sig must be positive!")
  }

  ## sequencing depth
  if (is.null(depth)){
    cat("Sequencing depths not given. Estimated by total counts.\n")
    depth <- colMeans(countmat)
  }
  depth <- exp(log(depth) - mean(log(depth)))

  ## start computing
  res <- matrix(NA, nrow=nrow(countmat), ncol=8)
  for (i in 1 : nrow(countmat))
  {
    res[i, ] <- robtest(y1=countmat[i, classlab == 1], d1=depth[classlab == 1],
                        y2=countmat[i, classlab == 2], d2=depth[classlab == 2],
                        c.tukey.beta=c.tukey.beta, c.tukey.phi=c.tukey.phi)
    if (i %% 10 == 0){
      cat(i, "genes processed.\n")
    }
  }

  ## convert the results into a data frame
  res <- data.frame(res)
  colnames(res) <- c('phi1', 'phi2', 'beta1', 'beta2', 'statistic.phi', 'statistic.beta',
                    'p.value.phi', 'p.value.beta')

  ## convert to fdr values
  res$fdr.phi <- p.adjust(res$p.value.phi, method="fdr")

  # FDR corrected p-values for DE test
  res$fdr.beta <- p.adjust(res$p.value.beta, method="fdr")

  return(res)
}
