#' \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function calculates \insertCite{jones1991;textual}{WData} kernel density estimator given a biased dataset and its bias function.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of `y`.
#' By default, it is set to the length-biased function.
#' @param x A numeric vector specifying the points where the density is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The smoothing bandwidth to be used in the density estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, options available are [`bw.f.BGMnrd0`][WData::bw.f.BGMnrd0()], [`bw.f.BGMcv`][WData::bw.f.BGMcv()],  [`bw.f.BGMboot1`][WData::bw.f.BGMboot1()] and [`bw.f.BGMboot2`][WData::bw.f.BGMboot2()].
#' Default is [`bw.f.BGMnrd0`][WData::bw.f.BGMnrd0()]. The specified (or computed) value of `bw` is multiplied by `adjust`.
#' @param adjust A numeric value representing the manual adjustment factor for the bandwidth. Default is 1.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound for  evaluation when `x` is not provided. Default is calculated based on the range of input data.
#' @param to Numeric value specifying the upper bound for  evaluation when `x` is not provided. Default is calculated based on the range of input data.
#' @param nb An integer specifying the number of points for evaluation when `x` is not provided. Default is 512.
#' @param plot A logical value indicating whether to plot the density estimation. Default is `TRUE`.
#' @param ... Additional arguments to be passed to bandwidth selection functions.
#' @return A list with the following components:
#'   \item{`x`}{The points where the density is estimated.}
#'   \item{`est_values`}{The estimated density values.}
#'   \item{`bw`}{The bandwidth used.}
#'   \item{`n`}{The sample size after elimination of missing values.}
#'   \item{`call`}{The call which produced the result.}
#'   \item{`data.name`}{The deparsed name of the y argument (biased dataset).}
#'   \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NA` values.}
#' @details
#' \insertCite{jones1991;textual}{WData} kernel density estimator is expressed as
#' \deqn{\widehat{f}_{\mathrm{J}}(y)=\frac{\widehat{\mu}_w}{n}\sum_{i=1}^{n}  \frac{1}{w(Y_i)} K_h(y-Y_i),\quad
#' \text{where}\quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n} \frac{1}{w(Y_i)}\right)^{-1},}
#' \eqn{h} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h}(u)=1/h K\left(u / h\right)}.
#' @references \insertAllCited{}
#' @seealso [`bw.f.BGMnrd0`][WData::bw.f.BGMnrd0()], [`bw.f.BGMcv`][WData::bw.f.BGMcv()],  [`bw.f.BGMboot1`][WData::bw.f.BGMboot1()] , [`bw.f.BGMboot2`][WData::bw.f.BGMboot2()]
#' @examples
#' # Rule of the thumb
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGMnrd0")
#' # Cross Validation
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGMcv")
#' # Bootstrap
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGMboot1", bw0 = "RT")
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGMboot1", bw0 = "Opt")
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGMboot2", nh = 50L)
df.jones <- function(y,
                     w = function(y) {
                       ifelse(y >= 0, y, NA)
                     },
                     x,
                     bw = "bw.f.BGMnrd0",
                     adjust = 1,
                     kernel = c(
                       "gaussian", "epanechnikov",
                       "rectangular", "triangular",
                       "biweight", "cosine", "optcosine"
                     ),
                     from,
                     to,
                     nb = 512L,
                     plot = TRUE,
                     ...) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  list2env(.get_xaxn_grid(y, x, from, to, nb, plot), envir = environment())

  # Bandwidth
  if (adjust <= 0) stop("'adjust' is not positive")

  if (is.character(bw)) {
    if (n < 2) stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(bw,
      bw.f.BGMnrd0 = bw.f.BGMnrd0(y, w, kernel),
      bw.f.BGMcv = bw.f.BGMcv(y, w, kernel, plot = FALSE, ...),
      bw.f.BGMboot1 = bw.f.BGMboot1(y, w, kernel, ...),
      bw.f.BGMboot2 = bw.f.BGMboot2(y, w, kernel, plot = FALSE, ...),
      stop("unknown bandwidth rule")
    )
  }
  if (!is.finite(bw)) stop("non-finite 'bw'")
  bw <- adjust * bw #  Manual Adjustment on the bw
  if (bw <= 0) stop("'bw' is not positive")

  # Density Estimation
  aux <- bw^(-1) * outer(x, y, "-")
  aux <- kernel_function_density(aux)
  aux <- aux %*% diag(weights)
  yords <- uw / bw * rowMeans(aux)

  if (plot == TRUE) {
    ord <- order(x)
    plot(x[ord], yords[ord],
      type = "l", main = "Jones' Estimator",
      xlab = paste(
        "n = ", n, "Bandwidth = ",
        format(round(bw, 5), nsmall = 5)
      ),
      ylab = "Density", col = "blue"
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
