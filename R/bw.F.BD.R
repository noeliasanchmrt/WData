#' \insertCite{bose2022;textual}{WData} bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function calculates the bandwidth for \insertCite{bose2022;textual}{WData} distribution estimator following the method proposed by \insertCite{bose2022;textual}{WData}.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function to be used in the estimation. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param x A numeric vector containing the points for which the bandwidth is calculated.
#' @param c_adj A numeric vector representing the constants to be used in the bandwidth calculation for each point of `x`.
#' @return A numeric vector containing the bandwidths for each point in `x`.
#' @details \insertCite{bose2022;textual}{WData} propose using local bandwidths the form
#' \deqn{\widehat{h}_{\mathrm{BD}, C(y)} (y) = \frac{C(y) \widehat{\sigma}}{(n w(y))^{1/3}},}
#' where \eqn{C(y)} is a positive function that depends on the point \eqn{y}, and \eqn{\widehat{\sigma}} is an estimate of the standard deviation of the distribution. \insertCite{borrajo2017;textual}{WData} suggests two options for estimating the standard deviation \eqn{\sigma}. The first is to take the square root of the following estimator of \eqn{\sigma^2}:
#' \deqn{\widehat{\sigma}_w^2=\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i
#' )}\right)^{-1}\right].}
#' The second option is to consider a robust estimator of \eqn{\sigma}, such as
#' \deqn{\widehat{\sigma}_{\mathrm{IQR}}=\mathrm{IQR}/(\Phi(0.75)-\Phi(0.25)),}
#' where \eqn{\mathrm{IQR}} is the interquartile range of the sample \eqn{\{Y_1, Y_2, \dots, Y_n\}} and \eqn{\Phi} is the cumulative distribution function of the standard normal distribution.
#' In order to prevent oversmoothing in the estimation, this
#' function considers \eqn{\widehat{\sigma}= \min\{\widehat{\sigma}_w, \widehat{\sigma}_{\mathrm{IQR}}\}}.
#' If some bandwidths are not positive, they are replaced by the mean of the neighbors.
#'
#' The simulations carried out by \insertCite{bose2022;textual}{WData} suggest that choosing \eqn{C(y)=0.25} or \eqn{C(y)=0.5} provides good results in the tail region of the distribution, with tails defined as points below the 5th percentile or above the 95th percentile. On the other hand, \eqn{C(y)=1.3} provides good results for the remaining points.
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @examples
#' bw.F.BD(shrub.data$Width, x = seq(0, 1, length.out = 512), c_adj = rep(1, 512))
bw.F.BD <- function(y,
                    w = function(y) {
                      ifelse(y >= 0, y, NA)
                    },
                    x,
                    c_adj) {
  list2env(.check_biased_dataset(y, w), envir = environment())

  if (!is.vector(x)) {
    stop("argument 'x' must be a vector")
  }

  if (missing(x)) {
    stop("argument 'x' is missing")
  }

  if (any(is.null(x))) {
    stop("argument 'x' must not contain NULL values")
  }

  if (all(c_adj <= 0)) {
    stop("argument 'c_adj' must be positive")
  }

  if (!is.vector(c_adj)) {
    stop("argument 'c_adj' must be a vector")
  }

  if (any(!is.numeric(c_adj))) {
    stop("argument 'c_adj' must be numeric")
  }

  sigma <- min(sqrt(uw * (mean(yw) - uw)), IQR(y) / 1.34)
  bw <- c_adj * sigma * n^(-0.33) * x^(-0.33)

  if (any(is.nan(bw)) || any(is.na(bw)) || any(bw <= 0)) {
    warning("some 'bw' are not positive. They will be replaced by the mean of the neighbors")
  }

  replace_with_closest_non_na <- function(bw) {
    bw_filled <- bw

    for (i in 1:length(bw)) {
      if (is.na(bw_filled[i]) || is.infinite(bw_filled[i]) || is.nan(bw_filled[i]) || bw_filled[i] <= 0) {
        left_value <- NA
        if (i > 1) {
          for (j in (i - 1):1) {
            if (is.finite(bw_filled[j])) {
              left_value <- bw_filled[j]
              break
            }
          }
        }

        right_value <- NA
        if (i < length(bw)) {
          for (j in (i + 1):length(bw)) {
            if (is.finite(bw_filled[j])) {
              right_value <- bw_filled[j]
              break
            }
          }
        }

        if (!is.na(left_value) && !is.na(right_value)) {
          bw_filled[i] <- mean(c(left_value, right_value), na.rm = TRUE)
        } else if (!is.na(left_value)) {
          bw_filled[i] <- left_value
        } else if (!is.na(right_value)) {
          bw_filled[i] <- right_value
        } else {
          bw_filled[i] <- bw_filled[i]
        }
      }
    }
    return(bw_filled)
  }

  bw <- replace_with_closest_non_na(bw)
}
