#' Rule of the thumb bandwidth selection for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function calculates the bandwidth for \insertCite{jones1991;textual}{WData} density estimator using the rule of the thumb.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The optimal bandwidth for \insertCite{jones1991;textual}{WData} density estimator based on the rule of the thumb.
#' @details The bandwidth \eqn{h} that minimizes the asymptotic integrated mean squared error when the kernel is density function \eqn{K} is given by
#' \deqn{h_{\mathrm{AMISE}}=\left(\frac{R(K) \mu_w \bar{\mu}_w}{n \sigma_K^4 R\left(f^{(2)}\right)}\right)^{1 / 5},}
#' where \eqn{R(K)} and \eqn{\sigma_K^4} depend only on the kernel function.
#' It would be assumed that \eqn{f} is normally distributed with mean \eqn{\mu} and standard deviation \eqn{\sigma} so that
#' \eqn{R\left(f^{(2)}\right)=3/8\pi^{-1/2}\sigma^{-5}}. \eqn{\mu_w} and \eqn{\bar{\mu}_{w}} are estimated as follows:
#' \deqn{\widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' \quad \text{and} \quad
#' \widehat{\bar{\mu}}_w= \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}.}
#' \insertCite{borrajo2017;textual}{WData} suggests two options for estimating \eqn{\sigma}.
#' The first one is to take
#' \deqn{\widehat{\sigma}_w=\sqrt{\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}
#' \left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\right]}.}
#' The second option is to consider a robust estimator of \eqn{\sigma} such as, for example,
#' \deqn{\widehat{\sigma}_{\mathrm{IQR}}=\mathrm{IQR}/(\Phi(0.75)-\Phi(0.25)),}
#' where IQR is the interquartile range and \eqn{\Phi} is the standard normal distribution function.
#' The bandwidth would then be obtained as
#' \deqn{\widehat{h}_{RT} = \left(\frac{8 \sqrt{\pi} R(K) \widehat{\mu}_w \widehat{\bar{\mu}}_w }
#' {3n \sigma_K^4}\right)^{1 / 5}\widehat{\sigma},}
#' where \eqn{\widehat{\sigma}} would be \eqn{\widehat{\sigma}_w^2} or \eqn{\widehat{\sigma}_{\mathrm{IQR}}}.
#' In order to prevent oversmoothing in the estimation, this
#' function considers \eqn{\widehat{\sigma}= \min\{\widehat{\sigma}_w, \widehat{\sigma}_{\mathrm{IQR}}\}}.
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @examples
#' bw.f.BGMnrd0(shrub.data$Width)
#' bw.f.BGMnrd0(shrub.data$Width, kernel = "epanechnikov")
bw.f.BGMnrd0 <- function(y,
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

  sigma <- sqrt(uw * (mean(yw) - uw))
  sigma * (8 * sqrt(pi) * RK * uw * uwb)^(0.2) * (3 * n * sigma_K_2^2)^(-0.2)
}
