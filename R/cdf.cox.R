#'  \insertCite{cox2005;textual}{WData} distribution estimator
#'
#' This function computes \insertCite{cox2005;textual}{WData} distribution estimator given a sample and the corresponding biased function.
#'
#' @param y  A numeric vector containing the biased sample.
#' @param w  A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @return A function of class `ecdf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details \insertCite{cox2005;textual}{WData} distribution estimator is expressed as
#' \deqn{\widehat{F}_n(y) = \frac{\widehat{\mu}_w}{n}\sum_{i=1}^{n} \frac{1}{w(Y_i)} \mathbb{I}(Y_i \leq y),
#' \quad \text{where} \quad
#' \widehat{\mu}_w = n \left(\sum_{i=1}^{n} \frac{1}{w(Y_i)} \right)^{-1}.}
#' @references \insertAllCited{}
#' @export
#' @examples
#' cdf.cox(y = shrub.data$Width)
cdf.cox <- function(y,
                    w = function(y) {
                      ifelse(y >= 0, y, NA)
                    }) {
  list2env(.check_biased_sample(y, w), envir = environment())

  vals <- unique(sort(y))
  weightsvals <- unique(sapply(vals, w)^(-1))
  rval <- approxfun(vals, uw * cumsum(weightsvals * tabulate(match(y, vals))) / n,
    method = "constant", yleft = 0, yright = 1,
    f = 0, ties = "ordered"
  )
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
