#' compute the empirical cumulative density.
#'
#' compute the empirical cumulative density given two random variables,
#' which will be the input for the mixture copula model
#'
#' @param x1 a vector of data for replicate 1, i.g., sequencing counts for each genomic locus.
#' @param x2 a vector of data for replicate 2. It must have the same length with x1.
#' @details The data was first ranked with a ties method equals to "random". The cumulative density
#' function (cdf) is then obtained using the ecdf function. To avoid infinity, a factor is multiplied
#' to the cdf. The factor is the \code{length(x1)/[length(x1)+1]}.
#' @return The empdist function returns a matrix that has two columns, each stores the empirical
#' cumulative density of one replicate.
#' @references  Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring reproducibility of high-throughput
#' experiments. Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
#' @export
#' @examples
#' data(Chipseq_TF)
#' x1 <- Chipseq_TF[,1]
#' x2 <- Chipseq_TF[,2]
#' empdist(x1, x2)

empdist = function(x1, x2){
  x1r = rank(x1, ties.method="random")
  x2r = rank(x2, ties.method="random")
  x1.cdf.func = ecdf(x1r); x2.cdf.func = ecdf(x2r)
  afactor <- length(x1r)/(length(x1r) + 1)
  x1.cdf <- x1.cdf.func(x1r) * afactor
  x2.cdf <- x2.cdf.func(x2r) * afactor
  u = cbind(x1.cdf, x2.cdf)
  return(u)
}
