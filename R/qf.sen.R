#' \insertCite{sen1984;textual}{WData} quantile estimator
#'
#' This function computes \insertCite{sen1984;textual}{WData} quantile estimator given a sample and the corresponding biased function.
#'
#' @param y  A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of the sample `y`.
#' By default, it is set to the length-biased function.
#' @return A function of class `eqf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details \insertCite{sen1984;textual}{WData} quantile estimator is expressed as
#' \deqn{ \widehat{F^{-1}_{\mathrm{S}}}(\tau) = Y_{(i)},
#' \quad \text{where} \quad i = \max \left\{ i: \sum_{j=1}^{i} \frac{1}{Y_{(j)}} \bigg/ \sum_{k=1}^{n} \frac{1}{Y_{(k)}} \leq \tau  \right\}.}
#' If the sample size \eqn{n} is small and \eqn{\tau} is very close to zero,
#' the inequality may not hold for any \eqn{i}, in which case \eqn{i} is taken as 1.
#' @references \insertAllCited{}
#' @export
#' @examples
#' qf.sen(y = shrub.data$Width)
qf.sen <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   }) {
  list2env(.check_biased_sample(y, w), envir = environment())

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw
  rval <- approxfun(ti, vals,
    method = "constant", rule = 2,
    f = 0, ties = "ordered"
  )
  class(rval) <- c("eqf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
