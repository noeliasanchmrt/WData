#' Rule of thumb bandwidth selector for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function computes the bandwidth selector for \insertCite{jones1991;textual}{WData} density estimator using the rule of thumb.
#'
#' @param y A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @param kernel A character string specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The rule of thumb bandwidth value.
#' @details The bandwidth is given by
#' \deqn{\widehat{h}_{f,\mathrm{RT}}
#' = \left(\frac{8 \sqrt{\pi} \widehat{\mu}_{w} \widehat{\bar{\mu}}_{w} R(K)}{3n \eta(K)^{2}}\right)^{1 / 5}\widehat{\sigma}_{w},}
#' where \eqn{R(K)} and \eqn{\eta(K)} depend only on the kernel and are defined as
#' \deqn{
#' R(K) = \int_{-\infty}^{+\infty} K(u)^2 du
#' \quad \text{and} \quad
#' \eta(K) = \int_{-\infty}^{+\infty} u^2 K(u) du.
#' }
#' The estimators \eqn{\widehat{\mu}_w} and \eqn{\widehat{\bar{\mu}}_w} are given by
#' \deqn{
#' \widehat{\mu}_w = n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' \quad \text{and} \quad
#' \widehat{\bar{\mu}}_w = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}.
#' }
#' \eqn{\widehat{\sigma}_{w}} is an estimation of the standard deviation of the distribution given by
#' \deqn{\widehat{\sigma}_w=\sqrt{\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}
#' \left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\right]}.}
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @export
#' @examples
#' bw.f.BGM.rt(shrub.data$Width)
#' bw.f.BGM.rt(shrub.data$Width, kernel = "epanechnikov")
bw.f.BGM.rt <- function(y,
                        w = function(y) {
                          ifelse(y >= 0, y, NA)
                        },
                        kernel = c(
                          "gaussian", "epanechnikov",
                          "rectangular", "triangular",
                          "biweight", "cosine", "optcosine"
                        )) {
  list2env(.check_biased_sample(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  sigma <- sqrt(uw * (mean(yw) - uw))
  sigma * (8 * sqrt(pi) * kernel_r * uw * uwb)^(0.2) * (3 * n * kernel_eta^2)^(-0.2)
}
