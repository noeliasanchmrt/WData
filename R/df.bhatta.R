#' \insertCite{bhattacharyya1988;textual}{WData} density estimator
#'
#' This function computes \insertCite{bhattacharyya1988;textual}{WData} density estimator given a sample and the corresponding biased function.
#'
#' @param y A numeric vector containing the biased sample.
#' @param w A function representing the bias function to be used. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @param plot Logical indicating whether to plot the estimated density. Default is `TRUE`.
#' @param ... Additional arguments to be passed to [`density`][stats::density()] function as, for instance, `kernel` or `bw`:
#' * `kernel` A character string giving the kernel to be used. This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian", and may be abbreviated to a unique prefix (single letter).
#' * `bw` The smoothing bandwidth to be used in the density estimation. `bw` can also be a character string giving a rule to choose the bandwidth. Options available can be checked in [`bw.nrd`][stats::bw.nrd()]. Default is `"nrd0"`.
#' @return A list with the following components:
#'   \item{`y.seq`}{The points where the density is estimated.}
#'   \item{`f.hat`}{The estimated density values.}
#'   \item{`bw`}{The bandwidth value.}
#'   \item{`n`}{The sample size after removal of `NaN`, `Na` and `Inf`.}
#'   \item{`call`}{The call which produced the result.}
#'   \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NaN`, `Na` or `Inf`.}
#' @details \insertCite{bhattacharyya1988;textual}{WData} density estimator is computed as follows:
#' \deqn{\widehat{f}_{\mathrm{B}, h_{g}}(y)= \widehat{\mu}_w w(y)^{-1} \widehat{g}_{h_{g}}(y),
#' \quad \text{where} \quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n} \frac{1}{w(Y_i)}\right)^{-1},}
#' and \eqn{\widehat{g}_{h_{g}}(y)} is the kernel density estimate of the given data `y` using [`density`][stats::density()] function with main arguments `bw` and `kernel`.
#' @references \insertAllCited{}
#' @examples
#' # Rule of thumb
#' df.bhatta(shrub.data$Width, bw = "nrd0")
#' # Cross Validation
#' df.bhatta(shrub.data$Width, bw = "ucv")
#' # Sheather & Jones
#' bhata_sj <- df.bhatta(shrub.data$Width, bw = "SJ-ste")
#' # Rectangular kernel
#' df.bhatta(shrub.data$Width, bw = "nrd0", kernel = "epanechnikov")
#'
df.bhatta <- function(y,
                      w = function(y) {
                        ifelse(y >= 0, y, NA)
                      },
                      plot = TRUE,
                      ...) {
  list2env(.check_biased_sample(y, w), envir = environment())

  g_hat <- do.call(density, list(x = y, ...))
  f_hat <- uw * g_hat$y / w(g_hat$x)

  if (plot == TRUE) {
    order <- order(g_hat$x)
    plot(g_hat$x[order], f_hat[order],
      type = "l",
      ylim = c(0, max(max(g_hat$y, na.rm = TRUE), max(f_hat, na.rm = TRUE))),
      ylab = "Density",
      xlab = paste(
        "n = ", n, "Bandwidth = ",
        format(round(g_hat$bw, 5), nsmall = 5)
      ),
      col = "blue",
      main = "Bhattacharyya's Estimator"
    )
    rug(y)
    lines(x = g_hat$x, y = g_hat$y, col = "black")
    legend("topright",
      legend = c(
        "Estimated Biased Density",
        "Estimated Unbiased Density"
      ),
      col = c("black", "blue"),
      lty = c(1, 1), cex = 0.8
    )
  }

  list(
    y.seq = g_hat$x,
    f.hat = f_hat,
    bw = g_hat$bw,
    data.name = data.name,
    n = n,
    has.na = has.na,
    call = match.call()
  )
}
