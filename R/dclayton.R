#' calculate the density of Clayton copula.
#'
#' compute the density of Clayton copula given the two individual variables.
#'
#' @param u a vector containing cumulative densities for one replicate data.
#' @param v a vector containing cumulative densities for the other replicate data.
#' @param beta the single paramter for the Clayton copula. It has a range of (0, \eqn{\infty}),
#' perfect dependence is archieved if \eqn{\beta} approximates \eqn{\infty}, while \eqn{\beta} approximates 0
#' implies no dependence.
#' @details The formula for the density of bivariate Clayton copula is:
#' \deqn{c(u, v | \beta) = (1+\beta)*((u*v)^(-1-\beta))*(u^(-\beta)+v^(-\beta)-1)^(-1/\beta-2)}
#' where \eqn{\beta > 0} in our application.
#' @references Nelsen, R. (2006). An Introduction to Copula, Second Edition, Springer.
#' G. G. Venter (2001). Tails of copulas. In Proceedings ASTIN Washington, USA, pages 68-113.
#' @return dclayton calculates the density of Clayton copula given the cumulative densities of
#' two random variables and the parameter of Clayton copula. The empirical cumulative densities
#' of the two random variables could be obtained using the function \code{empdist()}. Invalid arguments
#' will result in value \code{NaN}.
#' @export
#' @examples
#' data(Chipseq_TF)
#' x1 <- Chipseq_TF[,1]
#' x2 <- Chipseq_TF[,2]
#' U=empdist(x1, x2)
#' u <- U[,1]
#' v <- U[,2]
#' beta <- 2
#' dclayton(u, v, beta)

dclayton = function(u, v, beta){
  if (beta <=0 && !is.numeric(beta))
    stop("invalid argument : beta\n")
  else{
    m = u^(-beta)+v^(-beta)-1
    dclay=(1+beta)*((u*v)^(-1-beta))*(u^(-beta)+v^(-beta)-1)^(-1/beta-2)
  }
  return(dclay)
}



