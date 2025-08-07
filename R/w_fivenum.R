#' Tukey five-number summary
#'
#' This function calculates Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for biased data.
#'
#' @param x	 Numeric, maybe including `NA`s and `Inf`s.
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of `x`. By default, it is set to the length-biased function.
#' @param qmethod Character string specifying the quantile estimator to be used. Options are `"qf.SBC"` ([`qf.SBC`][WData::qf.SBC()]) and `"qf.sen"` ([`qf.sen`][WData::qf.sen()]). Default is `"qf.SBC"`.
#' @param na.rm Logical; if `TRUE`, all `NA`s and `NaN`s are dropped, before the statistics are computed.
#' @return A numeric vector of length 5 containing the Tukey five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for biased data.
#' @details This function is a modification of the [`fivenum`][stats::fivenum()] function, which is used to calculate Tukey's five number summary. The main difference is that this function allows to compute the five-number summary for biased data by calculating the five-number summary of the data based on [`qf.sen`][WData::qf.sen()] quantile estimator or in [`qf.SBC`][WData::qf.SBC()] quantile estimator.
#' @seealso [`w_boxstats`][WData::w_boxstats()].
#' @examples
#' w_fivenum(shrub.data$Width)
w_fivenum <- function(x, w = function(y) {
                        ifelse(y >= 0, y, NA)
                      },
                      qmethod = c("qf.SBC", "qf.sen"),
                      na.rm = TRUE) {
  list2env(.check_biased_dataset(x, w), envir = environment())

  xna <- is.na(x)
  if (any(xna)) {
    if (na.rm) {
      x <- x[!xna]
    } else {
      return(rep.int(NA, 5))
    }
  }
  x <- sort(x)
  n <- length(x)
  if (n == 0) {
    rep.int(NA, 5)
  } else {
    qmethod <- match.arg(qmethod)

    switch(qmethod,
      qf.SBC = c(min(x), qf.SBC(y = x, w = w)(c(0.25, 0.5, 0.75)), max(x)),
      qf.sen = c(min(x), qf.sen(y = x, w = w)(c(0.25, 0.5, 0.75)), max(x)),
      stop("unknown method for quantile computation")
    )
  }
}
