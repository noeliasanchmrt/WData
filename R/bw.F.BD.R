#' \insertCite{bose2022;textual}{WData} local bandwidth selector for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function implements the local bandwidth selector proposed by \insertCite{bose2022;textual}{WData} for their own kernel distribution estimator.
#'
#' @param y A numeric vector containing the biased sample.
#' @param w A function representing the bias function to be used. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @param y.seq A numeric vector containing the points on which the local bandwidth is estimated.
#' @param cy.seq A numeric vector representing the constants to be used in the bandwidth estimation for each point of `y.seq`.
#' Alternatively, a single numeric value can be provided, which will be used for all points in the `y.seq` vector.
#' @return A numeric vector containing the bandwidths for each point in `y.seq`.
#' @details Local bandwidths selectors are estimated using the formula:
#' \deqn{\widehat{h}_{F, \mathrm{BD}, C(y)} (y) = \frac{C(y) \widehat{\sigma}_{w}}{(n w(y))^{1/3}},}
#' where \eqn{C(y)} is a positive parameter that depends on the point \eqn{y} and \eqn{\widehat{\sigma}_{w}} is an estimation of the standard deviation of the distribution given by
#' \deqn{\widehat{\sigma}_w=\sqrt{\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}
#' \left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\right]}.}
#'
#' The parameter \eqn{C(y)} is provided to the function by using the argument `cy.seq`, which is a vector of positive values that is used to compute the bandwidth for each point in `y.seq`. Alternatively, a single numeric value can be provided, which will be used for all points in the `y.seq` vector.
#' The simulations carried out by \insertCite{bose2022;textual}{WData} suggest that choosing \eqn{C(y)=0.25} or \eqn{C(y)=0.5} provides good results in the tail region of the distribution, with tails defined as points below the 5th percentile or above the 95th percentile. On the other hand, \eqn{C(y)=1.3} provides good results for the remaining points.
#' If some bandwidths are not positive, they are replaced by the mean of the neighbors.
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @examples
#' bw.F.BD(shrub.data$Width, y.seq = seq(0, 1, length.out = 512), cy.seq = rep(1, 512))
bw.F.BD <- function(y,
                    w = function(y) {
                      ifelse(y >= 0, y, NA)
                    },
                    y.seq,
                    cy.seq) {
  list2env(.check_biased_sample(y, w), envir = environment())

  if (!is.vector(y.seq)) {
    stop("argument 'y.seq' must be a vector")
  }

  if (missing(y.seq)) {
    stop("argument 'y.seq' is missing")
  }

  if (any(is.null(y.seq))) {
    stop("argument 'y.seq' must not contain NULL values")
  }

  if (all(cy.seq <= 0)) {
    stop("argument 'cy.seq' must be positive")
  }

  if (!is.vector(cy.seq)) {
    stop("argument 'cy.seq' must be a vector")
  }

  if (any(!is.numeric(cy.seq))) {
    stop("argument 'cy.seq' must be numeric")
  }

  if (length(cy.seq) > 1 & length(cy.seq) != length(y.seq)) {
    stop("arguments 'y.seq' and 'cy.seq' must have the same length or 'cy.seq' must be a single value")
  }

  sigma <- sqrt(uw * (mean(yw) - uw))
  bw <- cy.seq * sigma * n^(-1 / 3) * w(y.seq)^(-1 / 3)


  if (any(is.nan(bw)) || any(is.na(bw)) || any(bw <= 0)) {
    warning("some 'bw' are not positive or could not be computated. They will be replaced by the mean of the neighbors")
  }

  # In case of NA, NaN or negative values, replace them with the mean of the closest non-NA, finite and positive values.
  # This is done by sorting the values based on y.seq and then replacing them in the original order.
  # The intent is to do not produce errors when bw are extremely small do to a small cy.seq, big n or big y.seq.

  replace_with_closest_non_na <- function(bw) {
    ord <- order(y.seq)
    bw_sorted <- bw[ord]

    for (i in 1:length(bw)) {
      if (is.na(bw_sorted[i]) || is.infinite(bw_sorted[i]) || is.nan(bw_sorted[i]) || bw_sorted[i] <= 0) {
        left_value <- NA
        if (i > 1) {
          for (j in (i - 1):1) {
            if (is.finite(bw_sorted[j])) {
              left_value <- bw_sorted[j]
              break
            }
          }
        }

        right_value <- NA
        if (i < length(bw)) {
          for (j in (i + 1):length(bw)) {
            if (is.finite(bw_sorted[j])) {
              right_value <- bw_sorted[j]
              break
            }
          }
        }

        if (!is.na(left_value) && !is.na(right_value)) {
          bw_sorted[i] <- mean(c(left_value, right_value), na.rm = TRUE)
        } else if (!is.na(left_value)) {
          bw_sorted[i] <- left_value
        } else if (!is.na(right_value)) {
          bw_sorted[i] <- right_value
        } else {
          bw_sorted[i] <- bw_sorted[i]
        }
      }
    }
    bw_filled <- numeric(length(bw))
    bw_filled[ord] <- bw_sorted
    return(bw_filled)
  }

  bw <- replace_with_closest_non_na(bw)
  return(bw)
}
