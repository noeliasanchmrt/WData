#' \insertCite{borrajo2017;textual}{WData} bootstrap bandwidth selection for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function calculates the bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator using the bias-corrected bootstrap method developed by \insertCite{borrajo2017;textual}{WData} and
#' based on the methodology introduced by \insertCite{cao1990;textual}{WData} and \insertCite{cao93;textual}{WData}.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param bw0 A character string specifying the method to determine the pilot bandwidth. Options are `"RT"` for rule of thumb and `"Opt"` for optimal bandwidth. Default is `"RT"`.
#' @return The estimated bootstrap bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator.
#' @details As stated by \insertCite{borrajo2017;textual}{WData}, optimal bootstrap bandwidth is given by
#' \deqn{
#' h_{\mathrm{AMISE}^*}= \left( \frac{R(K^{\ast})\widehat{\mu}_w\widehat{\bar{\mu}}_w}{n\sigma^4_{K^{\ast}} R\left(\widehat{f}_{J,h_0}^{(2)}\right)}\right)^{1/5},
#' \quad \text{where}  \quad
#' \widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1},
#' \quad
#' \widehat{\bar{\mu}}_w= \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n} \frac{1}{w(Y_i)^2}}
#' and \eqn{R(K^{\ast})} and \eqn{\sigma^4_{K^{\ast}}} depend only on the kernel function.
#' A simple proposal for pilot bandwidth \eqn{h_0} is to take the one provided by the rule of thumb multiplied by a factor of \eqn{\sqrt[5]{n}/\sqrt[7]{n}}
#' to avoid oversmoothing problems. We would obtain:
#' \deqn{
#' \widehat{h}_{B_{\mathrm{RT}}}= \left( \frac{R(K^{\ast})\widehat{\mu}_w\widehat{\bar{\mu}}_w}{n\sigma^4_{K^{\ast}} R\left(\widehat{f}_{J,\widehat{h}_{0,\mathrm{RT}}}^{(2)}\right)}\right)^{1/5},
#' \quad
#' \text{where} \quad
#' \widehat{h}_{0, \mathrm{RT}}= \frac{n^{1/5}}{n^{1/7}} \widehat{h}_{\mathrm{RT}}.
#' }
#' This is implemented when `bw0="RT"`.
#' An alternative is to consider the following asymptotic expression for the optimal pilot bandwidth:
#' \deqn{
#' \widehat{h}_{B_{opt}}= \left( \frac{R(K^{\ast})\widehat{\mu}_w\widehat{\bar{\mu}}_w}{n\sigma^4_{K^{\ast}} R\left(\widehat{f}_{J,\widehat{h}_{0,opt}}^{(2)}\right)}\right)^{1/5},
#' \quad
#' \text{where} \quad
#' \widehat{h}_{0, opt}= \left(\frac{5 \widehat{\mu}_w\widehat{\bar{\mu}}_w R\left(L^{(2)}\right)  }{2 n \sigma_{L}^{2}  R\left(f^{(3)}\right) } \right)^{1 / 7}
#' }
#' and \eqn{R(f^{(3)})} is estimated by assuming \eqn{f} is normal. This is implemented by stating `bw0="Opt"`.
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @examples
#' # Estimate bandwidth using bootstrap method with default pilot bandwidth
#' bw.f.BGMboot1(y = shrub.data$Width)
#' # Estimate bandwidth using bootstrap method with optimal pilot bandwidth
#' bw.f.BGMboot1(y = shrub.data$Width, bw0 = "Opt")
bw.f.BGMboot1 <- function(y,
                          w = function(y) {
                            ifelse(y >= 0, y, NA)
                          },
                          kernel = c(
                            "gaussian", "epanechnikov",
                            "rectangular", "triangular",
                            "biweight", "cosine", "optcosine"
                          ),
                          bw0 = c("RT", "Opt")) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  # Pilot bandwidth

  if (is.character(bw0)) {
    sigma <- sqrt(uw * (mean(yw) - uw))

    bw0 <- match.arg(bw0)

    bw0 <- switch(bw0,
      RT = {
        if (kernel == "rectangular") stop("rectangular kernel is not supported for automatic bandwidth selection with RT pilot bandwidth")
        if (kernel == "triangular") stop("triangular kernel is not supported for automatic bandwidth selection with RT pilot bandwidth")
        n^(1 / 5 - 1 / 7) * sigma * (8 * sqrt(pi) * RK * uw * uwb)^(0.2) * (3 * n * sigma_K_2^2)^(-0.2)
      },
      Opt = {
        if (kernel == "rectangular") stop("rectangular kernel is not supported for automatic bandwidth selection with optimal pilot bandwidth")
        if (kernel == "triangular") stop("triangular kernel is not supported for automatic bandwidth selection with optimal pilot bandwidth")
        Rfprime3 <- 15 / (16 * sqrt(pi) * sigma^7)
        ((5 * uw * uwb * RKprime2) / (2 * n * sigma_K_2 * Rfprime3))^(1 / 7)
      },
      stop("unkown pilot bandwidth")
    )
  }
  if (!is.finite(bw0)) stop("non-finite 'bw0'")
  if (bw0 <= 0) stop("'bw0' is not positive")

  message(sprintf("Pilot Bandwidth for Bootstrap: %f", bw0))

  # Bootstrap bandwidth

  fJ_bw0_2_hat_2 <- function(z) {
    aux <- (z - y) / bw0
    aux <- kernel_function_density_deriv2(aux) %*% diag(weights)
    aux <- outer(aux, aux, "*")
    aux <- (uw / (n * bw0^3))^2 * sum(aux)
  }

  R_fJ_bw0_2_hat_2 <-
    integrate(Vectorize(fJ_bw0_2_hat_2),
      lower = -Inf,
      upper = +Inf,
      subdivisions = 1000, rel.tol = .Machine$double.eps^.15
    )$value

  h_B_hat <- (RK * uw * uwb)^(0.2) * (n * sigma_K_2^2 * R_fJ_bw0_2_hat_2)^(-0.2)

  return(h_B_hat)
}
