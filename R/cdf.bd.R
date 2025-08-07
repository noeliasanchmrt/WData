#' \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function calculates \insertCite{bose2022;textual}{WData} kernel distribution estimator given a biased dataset and its bias function.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param x A numeric vector specifying the points where the distribution is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The smoothing bandwidth to be used in the distribution estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, options available are [`bw.F.SBCnrd0`][WData::bw.F.SBCnrd0()], [`bw.F.BD`][WData::bw.F.BD()], [`bw.F.SBCcv`][WData::bw.F.SBCcv()] and [`bw.F.SBCpi`][WData::bw.F.SBCpi()]. Default is [`bw.F.SBCnrd0`][WData::bw.F.SBCnrd0()]. The specified (or computed) value of `bw` multiplied by `adjust`.
#' @param adjust A numeric value representing the manual adjustment factor for the bandwidth. Default is 1.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound for  evaluation when `x` not provided. Default is calculated based on the range of input data.
#' @param to Numeric value specifying the upper bound for  evaluation when `x` not provided. Default is calculated based on the range of input data.
#' @param nb An integer specifying the number of points for evaluation when `x` not provided. Default is 512.
#' @param plot A logical value indicating whether to plot the estimation. Default is `TRUE`.
#' @param correction A character string specifying the boundary correction to be applied. Options are "none", "left", "right" and "both". Default is "none".
#' @param ... Additional arguments to be passed to bandwidth selection functions.
#' @return A list with the following components:
#'   \item{`x`}{The points where the distribution is estimated.}
#'   \item{`est_values`}{The estimated distribution values.}
#'   \item{`bw`}{The bandwidth used.}
#'   \item{`n`}{The sample size after elimination of missing values.}
#'   \item{`call`}{The call which produced the result.}
#'   \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NA` values.}
#' @details
#' \insertCite{bose2022;textual}{WData} kernel distribution estimator is expressed as
#' \deqn{
#' \widehat{F}_{h_{F}}(y) = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^n \frac{1}{w(Y_i)} W_{h_{F}}\left(y-Y_i\right),
#' \quad
#' \text{where}\quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n}\frac{1}{w(Y_i)}\right)^{-1},}
#' \eqn{h_{F}} is the bandwidth, \eqn{W} is the kernel distribution function and \eqn{W_{h_{F}}(u) = W(u/{h_{F}})}.
#' @references \insertAllCited{}
#' @seealso [`bw.F.SBCnrd0`][WData::bw.F.SBCnrd0()], [`bw.F.BD`][WData::bw.F.BD()], [`bw.F.SBCcv`][WData::bw.F.SBCcv()], [`bw.F.SBCpi`][WData::bw.F.SBCpi()]
#' @examples
#' cdf.bd(shrub.data$Width, kernel = "epanechnikov")
#' cdf.bd(shrub.data$Width, bw = "bw.F.SBCcv")
cdf.bd <- function(y,
                   w = function(y) {
                     ifelse(y >= 0, y, NA)
                   },
                   x,
                   bw = "bw.F.SBCnrd0",
                   adjust = 1,
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
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  list2env(.get_xaxn_grid(y, x, from, to, nb, plot), envir = environment())

  if (adjust <= 0) {
    stop("'adjust' is not positive")
  }

  # Bandwidth
  if (is.character(bw)) {
    if (n < 2) {
      stop("need at least 2 points to select a bandwidth automatically")
    }
    bw <- switch(bw,
      bw.F.SBCnrd0 = bw.F.SBCnrd0(y, w, kernel),
      bw.F.BD = bw.F.BD(y, w, x, ...),
      bw.F.SBCcv = bw.F.SBCcv(y, w, kernel, plot = FALSE, ...),
      bw.F.SBCpi = bw.F.SBCpi(y, w, kernel, ...),
      stop("unknown bandwidth rule")
    )
  }

  if (!all(is.numeric(bw))) {
    stop("non-finite 'bw'")
  }

  bw <- adjust * bw ## Manual Adjustment on the bw

  if (is.vector(bw) & length(bw) == 1) {
    if (bw <= 0) stop("'bw' is not positive")
  }

  if (is.vector(bw) & length(bw) > 1) {
    if (any(bw <= 0, na.rm = TRUE)) {
      stop("some 'bw' are not positive")
    }
    if (length(bw) != length(x)) {
      stop("number of bandwidths must be equal to the number of points for evaluation")
    }
    if (any(is.nan(bw) | is.infinite(bw) | is.na(bw))) {
      warning("some 'bw' are not finite")
    }
  }

  # Distribution Estimation
  aux <- diag(bw^(-1), length(x)) %*% outer(x, y, "-") # We need for local bandwidth
  aux <- kernel_function_distribution(aux)
  aux <- aux %*% diag(weights)
  yords <- uw * apply(aux, 1, mean)

  # Boundary correction

  correction <- match.arg(correction)

  if (correction != "none") {
    if (correction == "left") {
      yords <- (yords - yords[1]) / (1 - yords[1])
    }
    if (correction == "right") {
      yords <- yords / tail(yords, 1)
    }
    if (correction == "both") {
      yords <- (yords - yords[1]) / (tail(yords, 1) - yords[1])
    }
  }

  if (plot == TRUE) {
    ord <- order(x)
    plot(
      x[ord],
      yords[ord],
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
    x = x,
    est_values = yords,
    data.name = data.name,
    bw = bw,
    n = n,
    has.na = has.na,
    call = match.call()
  )
}
