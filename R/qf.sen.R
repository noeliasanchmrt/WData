#' \insertCite{sen1984;textual}{WData} quantile estimator
#'
#' This function calculates \insertCite{sen1984;textual}{WData} quantile estimator.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @return A function of class `eqf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details
#' The estimator proposed by \insertCite{sen1984;textual}{WData} is defined as:
#' \deqn{ \widehat{F^{-1}_{\mathrm{SEN}}}(\tau) = Y_{(i)},
#' \quad \text{where} \quad i = \max \left\{ i: \sum_{j=1}^{i} \frac{1}{Y_{(j)}} \bigg/ \sum_{i=j}^{n} \frac{1}{Y_{(j)}} \leq \tau  \right\}.}
#' If the sample size \eqn{n} is small and \eqn{\tau} is very close to zero,
#' the inequality may not hold for any \eqn{i}, in which case \eqn{i} is taken as 1.
#' The estimator can also be written as:
#' \deqn{ \widehat{F^{-1}_{\mathrm{SEN}}} (\tau) = \sum_{i=1}^{n} Y_{(i)} \mathbb{I} \left( T_i \leq \tau < T_{i+1} \right),
#' \quad \text{where} \quad
#' T_i = \sum_{j=1}^{i} \frac{1}{Y_{(j)}} \bigg/ \sum_{i=j}^{n} \frac{1}{Y_{(j)}}.}
#' It is understood that \eqn{T_{n+1} = 1}.
#' @references \insertAllCited{}
#' @examples
#' qf.sen(y = shrub.data$Width)
qf.sen <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   }) {
  list2env(.check_biased_dataset(y, w), envir = environment())

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
