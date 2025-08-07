#' Rule of the thumb bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function calculates the bandwidth for  \insertCite{bose2022;textual}{WData} kernel distribution estimator using the rule of the thumb.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The optimal bandwidth for \insertCite{bose2022;textual}{WData} kernel distribution estimator based on the rule of the thumb.
#' @details
#' The bandwidth is computed as follows:
#' \deqn{\widehat{h}_{F, \mathrm{RT}} =
#' \left(\frac{4 \sqrt{\pi} \widehat{\mu}_{w} \widehat{\bar{\mu}}_{w}\kappa(K)}{n \eta(K)^{2} }\right)^{1/3} \widehat{\sigma}, }
#' where both \eqn{\kappa(K)} and \eqn{\eta(K)} depend only on the kernel,
#' \deqn{\widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' \quad \text{and} \quad
#' \widehat{\bar{\mu}}_w= \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}.}
#' \eqn{\widehat{\sigma}} is an estimate of the standard deviation of the distribution given by
#' \deqn{\widehat{\sigma}_w=\sqrt{\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}
#' \left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\right]}.}
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @examples
#' bw.F.SBCnrd0(shrub.data$Width)
#' bw.F.SBCnrd0(shrub.data$Width, kernel = "epanechnikov")
#'
bw.F.SBCnrd0 <- function(y,
                         w = function(y) {
                           ifelse(y >= 0, y, NA)
                         },
                         kernel = c(
                           "gaussian", "epanechnikov",
                           "rectangular", "triangular",
                           "biweight", "cosine", "optcosine"
                         )) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  # Compute the bandwidth
  sigma <- sqrt(uw * (mean(yw) - uw))
  sigma * (4 * sqrt(pi) * uw * uwb * kernel_kappa)^(1 / 3) * (n * kernel_eta^2)^(-1 / 3)
}
