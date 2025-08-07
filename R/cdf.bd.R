#' \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function computes \insertCite{bose2022;textual}{WData} kernel distribution estimator given a sample and the corresponding biased function.
#'
#' @param y  A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @param y.seq A numeric vector specifying the points where the distribution is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The bandwidth to be used in the distribution estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, options available are [`bw.F.SBC.rt`][WData::bw.F.SBC.rt()], [`bw.F.BD`][WData::bw.F.BD()], [`bw.F.SBC.cv`][WData::bw.F.SBC.cv()] and [`bw.F.SBC.pi`][WData::bw.F.SBC.pi()]. Default is [`bw.F.SBC.rt`][WData::bw.F.SBC.rt()].
#' @param kernel A character string specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound of the grid where the estimator is computed when `y.seq` is not provided. Default is computed based on the range of input data.
#' @param to Numeric value specifying the upper bound of the grid where the estimator is computed when `y.seq` is not provided. Default is computed based on the range of input data.
#' @param nb An integer specifying the number of points at which the estimator is computed when `y.seq` is not provided. Default is 512.
#' @param plot A logical value indicating whether to plot the estimation. Default is `TRUE`.
#' @param correction A character string specifying the boundary correction to be applied. Options are "none", "left", "right" and "both". Default is "none".
#' @param ... Additional arguments to be passed to bandwidth selection functions.
#' @return A list with the following components:
#'   \item{`y.seq`}{The points where the distribution is estimated.}
#'   \item{`F.hat`}{The estimated distribution values.}
#'   \item{`bw`}{The bandwidth value.}
#'   \item{`n`}{The sample size after removal of `NaN`, `Na` and `Inf`.}
#'   \item{`call`}{The call which produced the result.}
#'   \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NaN`, `Na` or `Inf`.}
#' @details
#' \insertCite{bose2022;textual}{WData} kernel distribution estimator is expressed as
#' \deqn{
#' \widehat{F}_{h_{F}}(y) = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^n \frac{1}{w(Y_i)} W_{h_{F}}\left(y-Y_i\right),
#' \quad
#' \text{where}\quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n}\frac{1}{w(Y_i)}\right)^{-1},}
#' \eqn{h_{F}} is the bandwidth, \eqn{W} is the kernel distribution function and \eqn{W_{h_{F}}(u) = W(u/{h_{F}})}.
#' \insertCite{bose2022;textual}{WData} propose a truncation correction for variables with compact support \eqn{[a,b]} for the estimator as follows:
#' \deqn{
#' \widehat{F}_{[a,b], h_{F}}(y)=\left\{
#' \begin{array}{ll}
#' 0, & y<a \\
#' \dfrac{\widehat{F}_{h_{F}}(y)-\widehat{F}_{h_{F}}(a)}{\widehat{F}_{h_{F}}(b)-\widehat{F}_{h_{F}}(a)}, & a \leq y<b \\
#' 1, & y \geq b.
#' \end{array}
#' \right.
#' }
#' The truncation correction is also valid for variables supported on \eqn{[a, +\infty)} or \eqn{(-\infty, b]}, replacing \eqn{\widehat{F}_{h_{F}}(b)} by 1 or \eqn{\widehat{F}_{h_{F}}(a)} by 0, respectively, in the above expression.  This correction is implemented in the `correction` argument, which can take values "none", "left", "right" or "both". If "left", the estimator is corrected to 0 for values less than the minimum of `y.seq`; if "right", it is corrected to 1 for values greater than the maximum of `y.seq`; if "both", it applies both corrections simultaneously.
#' @references \insertAllCited{}
#' @seealso [`bw.F.SBC.rt`][WData::bw.F.SBC.rt()], [`bw.F.BD`][WData::bw.F.BD()], [`bw.F.SBC.cv`][WData::bw.F.SBC.cv()], [`bw.F.SBC.pi`][WData::bw.F.SBC.pi()]
#' @export
#' @examples
#' cdf.bd(shrub.data$Width, kernel = "epanechnikov")
#' cdf.bd(shrub.data$Width, bw = "bw.F.SBC.cv")
cdf.bd <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   },
                   y.seq,
                   bw = "bw.F.SBC.rt",
                   kernel = c(
                     "gaussian",
                     "epanechnikov",
                     "rectangular",
                     "triangular",
                     "biweight",
                     "cosine",
                     "optcosine"
                   ),
                   from,
                   to,
                   nb = 512L,
                   plot = TRUE,
                   correction = c("none", "left", "right", "both"),
                   ...) {
  list2env(.check_biased_sample(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  list2env(.get_xaxn_grid(y, y.seq, from, to, nb, plot), envir = environment())

  # Bandwidth
  if (is.character(bw)) {
    if (n < 2) {
      stop("need at least 2 points to select a bandwidth automatically")
    }
    bw <- switch(bw,
      bw.F.SBC.rt = bw.F.SBC.rt(y, w, kernel),
      bw.F.BD = bw.F.BD(y, w, y.seq, ...),
      bw.F.SBC.cv = bw.F.SBC.cv(y, w, kernel, plot = FALSE, ...),
      bw.F.SBC.pi = bw.F.SBC.pi(y, w, kernel, ...),
      stop("unknown bandwidth rule")
    )
  }

  if (!all(is.numeric(bw))) {
    stop("non-finite 'bw'")
  }

  if (is.vector(bw) & length(bw) == 1) {
    if (bw <= 0) stop("'bw' is not positive")
  }

  if (is.vector(bw) & length(bw) > 1) {
    if (any(bw <= 0, na.rm = TRUE)) {
      stop("some 'bw' are not positive")
    }
    if (length(bw) != length(y.seq)) {
      stop("number of bandwidths must be equal to the number of points at which the estimator is computed")
    }
    if (any(is.nan(bw) | is.infinite(bw) | is.na(bw))) {
      warning("some 'bw' are not finite")
    }
  }

  # Distribution Estimation
  aux <- diag(bw^(-1), length(y.seq)) %*% outer(y.seq, y, "-") # We need for local bandwidth
  aux <- kernel_function_distribution(aux)
  aux <- aux %*% diag(weights)
  F.hat <- uw * apply(aux, 1, mean)

  # Boundary correction

  correction <- match.arg(correction)

  if (correction != "none") {
    if (correction == "left") {
      F.hat <- (F.hat - F.hat[1]) / (1 - F.hat[1])
    }
    if (correction == "right") {
      F.hat <- F.hat / tail(F.hat, 1)
    }
    if (correction == "both") {
      F.hat <- (F.hat - F.hat[1]) / (tail(F.hat, 1) - F.hat[1])
    }
  }

  if (plot == TRUE) {
    ord <- order(y.seq)
    plot(
      y.seq[ord],
      F.hat[ord],
      type = "l",
      main = "Bose & Dutta's Estimator",
      xlab = paste(
        "n = ", n, "Bandwidth = ",
        ifelse(is.vector(bw) && length(bw) > 1,
          "Local Bandwidth",
          format(round(bw, 5), nsmall = 5)
        )
      ),
      ylab = "Distribution",
      ylim = c(0, 1.1),
      col = "blue"
    )
    rug(y)
  }

  list(
    y.seq = y.seq,
    F.hat = F.hat,
    data.name = data.name,
    bw = bw,
    n = n,
    has.na = has.na,
    call = match.call()
  )
}
