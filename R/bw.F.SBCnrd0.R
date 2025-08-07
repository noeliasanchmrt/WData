#' Rule of the thumb bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function calculates the bandwidth for  \insertCite{bose2022;textual}{WData} kernel distribution estimator using the rule of the thumb.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The optimal bandwidth for \insertCite{bose2022;textual}{WData} kernel distribution estimator based on the rule of the thumb.
#' @details
#' The bandwidth that minimizes the Asymptotic Mean Integrated Squared Error (AMISE) is:
#' \deqn{h_{\mathrm{AMISE}} = \arg \min_{h>0} \mathrm{AMISE} (\widehat{F}_h)
#' = \left(\frac{\mu_{w} \bar{\mu}_{w} \tau_{W}^{2} }{n \sigma_K^4 R\left(f^{(1)}\right)}\right)^{1/3}, }
#' where \eqn{\tau_{W}^{2}} depends only on the kernel and \eqn{\mu_w} and \eqn{\bar{\mu}_{w}} are estimated as follows:
#' \deqn{\widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' \quad \text{and} \quad
#' \widehat{\bar{\mu}}_w= \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}.}
#' \insertCite{borrajo2017;textual}{WData} suggests two options for estimating \eqn{\sigma}.
#' The first one is to take
#' \deqn{\widehat{\sigma}_w=\sqrt{\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}
#' \left[\left(\frac{1}{n} \sum_{i=1}^n w(Y_i)\right)-\left(\frac{1}{n} \sum_{i=1}^n \frac{1}{w(Y_i)}\right)^{-1}\right]}.}
#' The second option is to consider a robust estimator of \eqn{\sigma} such as, for example,
#' \deqn{\widehat{\sigma}_{\mathrm{IQR}}=\mathrm{IQR}/(\Phi(0.75)-\Phi(0.25)),}
#' where IQR is the interquartile range and \eqn{\Phi} is the standard normal distribution function. If \eqn{f} is assumed to follow a normal distribution, then \eqn{R\left(f^{(1)}\right) = 1/(4 \sqrt{\pi} \sigma^3)},
#' where \eqn{\sigma} represents the standard deviation of the distribution and can be estimated by
#' \eqn{\widehat{\sigma}_w} or \eqn{\widehat{\sigma}_{\mathrm{IQR}}}. The smoothing parameter estimate is then given by:
#' \deqn{\widehat{h}_{\mathrm{RT}} = \left(\frac{4 \sqrt{\pi} \widehat{\mu}_w \widehat{\bar{\mu}}_{w} \tau_{W}^{2}}{n\sigma_{K}^{4}}\right)^{1/3} \widehat{\sigma},}
#' where \eqn{\widehat{\sigma}} represents \eqn{\widehat{\sigma}_w} or \eqn{\widehat{\sigma}_{\mathrm{IQR}}}.
#' In order to prevent oversmoothing in the estimation, this
#' function considers \eqn{\widehat{\sigma}= \min\{\widehat{\sigma}_w, \widehat{\sigma}_{\mathrm{IQR}}\}}.
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
  sigma * (4 * sqrt(pi) * uw * uwb * intudW2)^(1 / 3) * (n * sigma_K_2^2)^(-1 / 3)
}
