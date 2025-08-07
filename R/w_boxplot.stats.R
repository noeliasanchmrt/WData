#' Box plot statistics
#'
#' This function is typically called by another function to gather the statistics necessary for producing box plots, but may be invoked separately.
#'
#' @param x A numeric vector for which the boxplot will be constructed (`NA`s and `NaN`s are allowed and omitted).
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of `y`.
#' By default, it is set to the length-biased function.
#' @param qmethod Character string specifying the quantile estimator to be used. Options are `"qf.SBC"` ([`qf.SBC`][WData::qf.SBC()]) and `"qf.sen"` ([`qf.sen`][WData::qf.sen()]). Default is `"qf.SBC"`.
#' @param coef This determines how far the plot ‘whiskers’ extend out from the box. If `coef` is positive, the whiskers extend to the most extreme data point which is no more than `coef` times the length of the box away from the box. A value of zero causes the whiskers to extend to the data extremes (and no outliers be returned).
#' @param do.conf Logical; if `FALSE`, it will be empty in the result.
#' @param do.out Logical; if `FALSE`, it will be empty in the result.
#' @return A list with components:
#' \item{stats}{A vector of length 5, containing the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker.}
#' \item{n}{The number of non-`NA` observations in the sample.}
#' \item{conf}{The lower and upper extremes of the ‘notch’ (`if(do.conf)`). See the details on [`boxplot.stats`][WData::boxplot.stats()].}
#' \item{out}{The values of any data points which lie beyond the extremes of the whiskers (`if(do.out)`).}
#' @details This function is a modification of the [`boxplot.stats`][WData::boxplot.stats()] function, which is used to calculate the statistics needed to create a boxplot. The main difference is that this function allows to compute a boxplot for biased data by calculating the five-number summary of the data based on [`qf.sen`][WData::qf.sen()] quantile estimator or in [`qf.SBC`][WData::qf.SBC()] quantile estimator.
#' @seealso [`w_fivenum`][WData::w_fivenum()], [`w_boxplot`][WData::w_boxplot()].
#' @examples
#' w_boxstats(shrub.data$Width)
w_boxstats <- function(x, w = function(y) {
                         ifelse(y >= 0, y, NA)
                       },
                       qmethod = c("qf.SBC", "qf.sen"),
                       coef = 1.5, do.conf = TRUE, do.out = TRUE) {
  list2env(.check_biased_dataset(x, w), envir = environment())

  if (coef < 0) {
    stop("'coef' must not be negative")
  }
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- w_fivenum(x = x, w = w, qmethod = qmethod, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0) {
    do.out <- FALSE
  } else {
    out <- if (!is.na(iqr)) {
      x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * iqr)
    } else {
      !is.finite(x)
    }
    if (any(out[nna], na.rm = TRUE)) {
      stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
  }
  conf <- if (do.conf) {
    stats[3L] + c(-1.58, 1.58) * iqr / sqrt(n)
  }
  list(stats = stats, n = n, conf = conf, out = if (do.out) x[out & nna] else numeric())
}
