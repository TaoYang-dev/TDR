#' estimate the paramters of the mixture copula model.
#'
#' function to estimate the paramters of the mixture copula model (i.e., mixture of Gumbel and Clayton)
#' given the initialization parameters, maximum iteration and the stopping criterion.
#'
#' @param x1 a vector containing numeric values (i.g., sequencing counts of genomic loci) for one replicate data.
#' @param x2 a vector containing numeric values (i.g., sequencing counts of genomic loci) for the other replicate data.
#' @param p initialization of proportion paramter, any value between 0 and 1 (not include).
#' @param alpha initialization of Gumbel copula paramter, any value between 0 and 1 (not include).
#' @param beta  initialization of Clayton copula paramter, any value greater than 0 (not include).
#' @param it the maximum iterations. Usually 100 is enough to let estimation converge.
#' @param eps the stopping criterion. The default is 0.001, meaning if any of the paramter changes less than 0.001,
#' the estimation algorithm will stop.
#' @details The expectation maximization (EM) algorithm is adopted to estimate the paramters for the mixture copula model. The mixture copula model is:
#' \deqn{C(u, v | \pi, \alpha, \beta) = \pi * C_{clayton}(u, v | \beta) + (1-\pi) * C_{gumbel}(u, v | \alpha)}
#' where pi is the proportion paramters indicating how much of the mixture copula is contributed from Clayton copula, alpha is the single paramter for Gumbel copula, and the beta is the single paramter for Clayton copula.
#' @return
#'  \itemize{
#'   \item{trace.para} {the trace of the paramters for each iteration.}
#'   \item{trace.loglik} {the trace of the log likelihood for each iteration.}
#'   \item{membership} {a vector of binary data inidicating which data points are mapped to Clayton part (1) and
#'                      which are mapped to Gumbel copula (0).}
#'  }
#' @references Evaluating the Reproducibility and Quality of High Throughput Sequencing Data with Tail Dependences
#' of Mixture Copula (2016). Tao Yang, Feng Yue, Qunhua Li.
#' @export
#' @examples
#' data(Chipseq_TF)
#' x1 <- Chipseq_TF[,1]
#' x2 <- Chipseq_TF[,2]
#' TDR.est(x1, x2)

TDR.est = function(x1, x2, p=0.25, alpha=1.5, beta=1.5, it=100, eps=0.001) {

  start.time=proc.time()[3]

  U = empdist(x1, x2)
  u=U[,1];v=U[,2]

  if (p <=0 && p>=1 && !is.numeric(p))
    stop("invalid argument: intialize p with numeric number between 0 and 1 (not include)\n")

  else{

    p.trace = array()
    alpha.trace = array()
    beta.trace = array()

    loglik.gum = array()
    loglik.cla = array()
    loglik.T = array()

    i=0
    repeat {
      i=i+1

      d.gumb = dgumbel(u, v, alpha)
      d.clay = dclayton(u, v, beta)

      ki.p = p*d.clay/((1-p)*d.gumb+p*d.clay)
      p.p = sum(ki.p)/length(ki.p)

      fgum = function(alpha){
        sum((1-ki.p)*(log(1-p.p)+log(dgumbel(u, v, alpha))))
      }
      optima.gum = optimize(f = fgum, c(1 + sqrt(.Machine$double.eps), 10), maximum = T)
      alpha.p = optima.gum$maximum
      loglik.gum[i] = optima.gum$objective


      fclay = function(beta){
        sum(ki.p*(log(p.p)+log(dclayton(u, v, beta))))
      }
      optima.cla = optimize(f = fclay, c(0 + sqrt(.Machine$double.eps), 20), maximum = T)
      beta.p = optima.cla$maximum
      loglik.cla[i] = optima.cla$objective

      dif.p = abs(p.p-p)
      dif.al = abs(alpha.p-alpha)
      dif.be = abs(beta.p-beta)

      if(dif.p > eps){p = p.p} else {p = p}
      if(dif.al > eps){alpha = alpha.p} else {alpha = alpha}
      if(dif.be > eps){beta = beta.p} else {beta = beta}

      p.trace[i] = p
      alpha.trace[i] = alpha
      beta.trace[i] = beta

      loglik.T[i] = loglik.gum[i] + loglik.cla[i]

      if (dif.p < eps & dif.al < eps & dif.be < eps | i == it){break}
    }
    membership = ifelse(ki.p > 0.5, 1, 0)
    parameters = cbind(p.trace, alpha.trace, beta.trace)

    likelihood = cbind(loglik.gum, loglik.cla, loglik.T)
    est=list(parameters, likelihood, membership)
    names(est)=c("trace.para", "trace.loglik", "memebership")
  }
  print(proc.time()[3]-start.time)
  return(est)
}
