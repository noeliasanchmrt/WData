#' Plug-in bandwidth selector for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function computes the bandwidth selector for \insertCite{bose2022;textual}{WData} kernel distribution estimator using the plug-in method.
#'
#' @param y A numeric vector containing the biased sample.
#' @param w A function representing the bias function to be used. It must be evaluable and positive in each point of the sample the sample `y`. By default, it is set to the length-biased function.
#' @param kernel A character string specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The bandwidth value using the plug-in method.
#' @details The bandwidth is given by:
#' \deqn{\widehat{h}_{F, \mathrm{PI}} = \left(\frac{\widehat{\mu}_w \widehat{\bar{\mu}}_w \kappa(K) }{n \eta(K)^2 R \left(\widehat{f}_{\mathrm{J},\widehat{h}_{F,0, \mathrm{opt}}}^{(1)}\right)}\right)^{1/3},}
#' where \eqn{\kappa(K)} and \eqn{\eta(K)} depend only on the kernel and are defined as
#' \deqn{
#' \kappa(K) = \int_{-\infty}^{+\infty} 2 u W (u) K(u) du
#' \quad \text{and} \quad
#' \eta(K) = \int_{-\infty}^{+\infty} u^2 K(u) du,
#' }
#' where \eqn{ W} is the kernel distribution function associated with the kernel density function \eqn{K}.
#' The estimators \eqn{\widehat{\mu}_w} and \eqn{\widehat{\bar{\mu}}_w} are given by
#' \deqn{
#' \widehat{\mu}_w = n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' \quad \text{and} \quad
#' \widehat{\bar{\mu}}_w = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}.
#' }
#' \eqn{\widehat{h}_{F,0, \mathrm{opt}}} is an estimator of
#' \deqn{h_{F,0, \mathrm{opt}} = \left(\frac{3 \mu_w \bar{\mu}_w R\left(L^{(1)}\right)}{2 n \eta(L) R\left( f^{(2)} \right)^2}\right)^{1/5},}
#' where \eqn{R\left(f^{(2)}\right)} is estimated assuming that \eqn{f} follows a gaussian distribution and \eqn{\mu_w} and \eqn{\bar{\mu}_w} are estimated by \eqn{\widehat{\mu}_w} and \eqn{\widehat{\bar{\mu}}_w} as defined above.
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @export
#' @examples
#' bw.F.SBC.pi(shrub.data$Width, kernel = "epanechnikov")
bw.F.SBC.pi <- function(y,
                        w = function(y) {
                          ifelse(y >= 0, y, NA)
                        },
                        kernel = c(
                          "gaussian",
                          "epanechnikov",
                          "rectangular",
                          "triangular",
                          "biweight",
                          "cosine",
                          "optcosine"
                        )) {
  list2env(.check_biased_sample(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  # Pilot bandwidth
  if (kernel == "rectangular") {
    stop("rectangular kernel is not supported for plug-in pilot bandwidth")
  }
  sigma <- sqrt(uw * (mean(yw) - uw))
  bw0 <- sigma * (4 * sqrt(pi) * kernel_r_deriv1 * uw * uwb)^(0.2) * (n * kernel_eta)^(-0.2)
  if (!is.finite(bw0)) stop("non-finite 'bw0'")
  if (bw0 <= 0) stop("'bw0' is not positive")

  message(sprintf("Pilot Bandwidth: %f", bw0))

  # Bootstrap bandwidth

  fJ_bw0_1_hat_2 <- function(z) { # Squared first derivative of Jones' density estimator
    aux <- (z - y) / bw0
    aux <- kernel_function_density_deriv1(aux) %*% diag(weights)
    aux <- (uw / (n * bw0^2)) * sum(aux)
    aux <- aux^2
  }

  R_fJ_bw0_1_hat_2 <- integrate(Vectorize(fJ_bw0_1_hat_2),
    lower = -Inf,
    upper = +Inf,
    subdivisions = 1000, rel.tol = .Machine$double.eps^.15
  )$value

  ((kernel_kappa * uw * uwb) / (n * kernel_eta^2 * R_fJ_bw0_1_hat_2))^(1 / 3)
}
