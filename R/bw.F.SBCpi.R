#' Plug-in bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function calculates the bandwidth for \insertCite{bose2022;textual}{WData} kernel distribution estimator using the plug-in method.
#'
#' @param y A numeric vector for which the density estimation is performed.
#' @param w A function representing the bias function to be used in the estimation. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return The optimal bandwidth for \insertCite{bose2022;textual}{WData} kernel distribution estimator using the plug-in method.
#' @details The bandwidth is given by:
#' \deqn{\widehat{h}_{\mathrm{PI}} = \left(\frac{\widehat{\mu}_w \widehat{\bar{\mu}}_w \tau_{W}^{2} }{n \sigma_K^4 R \left(\widehat{f}_{\mathrm{J},\widehat{g}}^{(1)}\right)}\right)^{1/3},
#' \quad \text{where} \quad
#' \widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1},
#' \quad
#' \widehat{\bar{\mu}}_w= \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}}
#' and \eqn{\tau_{W}^{2}} depends only on the kernel function. \eqn{\widehat{g}} is an estimator of:
#' \deqn{g = \arg \min_{g>0} \mathrm{AMSE} \left(R\left(\widehat{f}_{\mathrm{J},g}^{(1)}\right)\right)
#' = \left(\frac{3 \mu_w \bar{\mu}_w R\left(L^{(1)}\right)}{2 n \sigma_L^2 R\left( f^{(2)} \right)^2}\right)^{1/5}.}
#' Here, \eqn{R\left(f^{(2)}\right)} is estimated assuming that \eqn{f} follows a normal distribution.
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @examples
#' bw.F.SBCpi(shrub.data$Width, kernel = "epanechnikov")
bw.F.SBCpi <- function(y,
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
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  # Pilot bandwidth
  if (kernel == "rectangular") {
    stop("rectangular kernel is not supported for optimal pilot bandwidth")
  }
  sigma <- sqrt(uw * (mean(yw) - uw))
  bw0 <- sigma * (4 * sqrt(pi) * RKprime * uw * uwb)^(0.2) * (n * sigma_K_2)^(-0.2)
  if (!is.finite(bw0)) stop("non-finite 'bw0'")
  if (bw0 <= 0) stop("'bw0' is not positive")

  message(sprintf("Pilot Bandwidth: %f", bw0))

  # Bootstrap bandwidth

  fJ_bw0_1_hat_2 <- function(z) { # Squared first derivative of Jones' density estimator
    aux <- (z - y) / bw0
    aux <- kernel_function_density_deriv(aux) %*% diag(weights)
    aux <- (uw / (n * bw0^2)) * sum(aux)
    aux <- aux^2
  }

  R_fJ_bw0_1_hat_2 <- integrate(Vectorize(fJ_bw0_1_hat_2),
    lower = -Inf,
    upper = +Inf,
    subdivisions = 1000, rel.tol = .Machine$double.eps^.15
  )$value

  ((intudW2 * uw * uwb) / (n * sigma_K_2^2 * R_fJ_bw0_1_hat_2))^(1 / 3)
}
