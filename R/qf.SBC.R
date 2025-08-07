#' Empirical quantile function based on \insertCite{cox2005;textual}{WData} distribution estimator
#'
#' This function implements the empirical quantile function based on the estimator
#' proposed by \insertCite{cox2005;textual}{WData} for the distribution function.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @return A function of class `eqf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details
#' The estimator is defined as:
#' \deqn{ \widehat{F^{-1}_n}(\tau) = Y_{(i)}
#' \quad \text{where} \quad
#' i = \min \left\{ i: \tau \leq  \sum_{j=1}^{i} \frac{1}{Y_{(j)}}
#' \bigg/ \sum_{i=1}^{n} \frac{1}{Y_{(j)}} \right\}.}
#' The estimator can also be written as:
#' \deqn{ \widehat{F^{-1}_n} (\tau) = \sum_{i=1}^{n} Y_{(i)} \mathbb{I} \left( T_{i-1} < \tau \leq T_{i} \right),
#' \quad \text{where}
#' \quad T_i = \sum_{j=1}^{i} \frac{1}{Y_{(j)}} \bigg/ \sum_{i=1}^{n} \frac{1}{Y_{(j)}}.}
#' It is understood that \eqn{T_{0} = 0}.
#' @references \insertAllCited{}
#' @seealso [`cdf.cox`][WData::cdf.cox()]
#' @examples
#' qf.SBC(y = shrub.data$Width)
#'
qf.SBC <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   }) {
  list2env(.check_biased_dataset(y, w), envir = environment())

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw
  rval <- approxfun(ti, vals,
    method = "constant", rule = 2,
    f = 1, ties = "ordered"
  )
  class(rval) <- c("eqf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
