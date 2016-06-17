#' extracting proportion, upper tail dependence and lower tail dependence
#'
#' function to extract proportion, upper tail dependence and lower tail dependence
#'
#' @param object the output from the function TDR.est
#' @return
#' \itemize{
#'  \item {pi} {the proportion parameter that indicates how much of mixture model is contributed by Clayton copula.}
#'  \item {UT} {the upper tail dependence coefficient.}
#'  \item {LT} {the lower tail dependence coefficient.}
#' }
#' @references Evaluating the Reproducibility and Quality of High Throughput Sequencing Data with Tail Dependences
#' of Mixture Copula (2016). Tao Yang, Feng Yue, Qunhua Li.
#' @export
#' @examples
#' data(Chipseq_TF)
#' x1 <- Chipseq_TF[,1]
#' x2 <- Chipseq_TF[,2]
#' U=empdist(x1, x2)
#' u <- U[,1]
#' v <- U[,2]
#' test <- TDR.est(u, v)
#' TailD(test)

TailD = function(object)
{
  trace=object[[1]]
  p= trace[nrow(trace),1]
  UT = (1-trace[nrow(trace),1])*(2-2^(1/trace[nrow(trace), 2]))
  LT = trace[nrow(trace),1]*2^(-1/trace[nrow(trace), 3])
  para = c(p, UT, LT)
  names(para) = c("pi", "UT", "LT")
  para
}
