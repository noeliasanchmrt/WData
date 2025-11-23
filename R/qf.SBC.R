#' Empirical quantile function based on \insertCite{cox2005;textual}{WData} distribution estimator
#'
#' This function computes the empirical quantile function based on the distribution function estimator
#' proposed by \insertCite{cox2005;textual}{WData}.
#'
#' @param y  A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of the sample `y`.
#' By default, it is set to the length-biased function.
#' @return A function of class `eqf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details
#' The estimator is defined as:
#' \deqn{ \widehat{F^{-1}_n}(\tau) = Y_{(i)}
#' \quad \text{where} \quad
#' i = \min \left\{ i: \tau \leq  \sum_{j=1}^{i} \frac{1}{Y_{(j)}}
#' \bigg/ \sum_{k=1}^{n} \frac{1}{Y_{(k)}} \right\}.}
#' @references \insertAllCited{}
#' @seealso [`cdf.cox`][WData::cdf.cox()]
#' @export
#' @examples
#' qf.SBC(y = shrub.data$Width)
#'
qf.SBC <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   }) {
  list2env(.check_biased_sample(y, w), envir = environment())

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
