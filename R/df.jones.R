#' \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function computes \insertCite{jones1991;textual}{WData} kernel density estimator given a sample and the corresponding biased function.
#'
#' @param y  A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points.
#' It must be evaluable and positive in each point of the sample `y`.
#' By default, it is set to the length-biased function.
#' @param y.seq A numeric vector specifying the points where the density is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The smoothing bandwidth to be used in the density estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, options available are [`bw.f.BGM.rt`][WData::bw.f.BGM.rt()], [`bw.f.BGM.cv`][WData::bw.f.BGM.cv()],  [`bw.f.BGM.boot1`][WData::bw.f.BGM.boot1()] and [`bw.f.BGM.boot2`][WData::bw.f.BGM.boot2()].
#' Default is [`bw.f.BGM.rt`][WData::bw.f.BGM.rt()].
#' @param kernel A character string specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound of the grid where the estimator is computed when `y.seq` is not provided. Default is computed based on the range of input data.
#' @param to Numeric value specifying the upper bound of the grid where the estimator is computed when `y.seq` is not provided. Default is computed based on the range of input data.
#' @param nb An integer specifying the number of points at which the estimator is computed when `y.seq` is not provided. Default is 512.
#' @param plot A logical value indicating whether to plot the density estimation. Default is `TRUE`.
#' @param ... Additional arguments to be passed to bandwidth selection functions.
#' @return A list with the following components:
#'   \item{`y.seq`}{The points where the density is estimated.}
#'   \item{`f.hat`}{The estimated density values.}
#'   \item{`bw`}{The bandwidth value.}
#'   \item{`n`}{The sample size after removal of `NaN`, `Na` and `Inf`.}
#'   \item{`call`}{The call which produced the result.}
#'   \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NaN`, `Na` or `Inf`.}
#' @details
#' \insertCite{jones1991;textual}{WData} kernel density estimator is expressed as
#' \deqn{\widehat{f}_{\mathrm{J}, h_{f}}(y)=\frac{\widehat{\mu}_w}{n}\sum_{i=1}^{n}  \frac{1}{w(Y_i)} K_{h_{f}}(y-Y_i),\quad
#' \text{where}\quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n} \frac{1}{w(Y_i)}\right)^{-1},}
#' \eqn{h_{f}} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h_{f}}(u)=1/h_{f} K\left(u / h_{f}\right)}.
#' @references \insertAllCited{}
#' @seealso [`bw.f.BGM.rt`][WData::bw.f.BGM.rt()], [`bw.f.BGM.cv`][WData::bw.f.BGM.cv()],  [`bw.f.BGM.boot1`][WData::bw.f.BGM.boot1()] , [`bw.f.BGM.boot2`][WData::bw.f.BGM.boot2()]
#' @examples
#' # Rule of thumb
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.rt")
#' # Cross Validation
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.cv")
#' # Bootstrap
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot1", bw0 = "RT")
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot1", bw0 = "PI")
#' df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot2", nh = 50L)
df.jones <- function(y,
                     w = function(y) {
                       ifelse(y >= 0, y, NA)
                     },
                     y.seq,
                     bw = "bw.f.BGM.rt",
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
  list2env(.check_biased_sample(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  list2env(.get_xaxn_grid(y, y.seq, from, to, nb, plot), envir = environment())

  if (is.character(bw)) {
    if (n < 2) stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(bw,
      bw.f.BGM.rt = bw.f.BGM.rt(y, w, kernel),
      bw.f.BGM.cv = bw.f.BGM.cv(y, w, kernel, plot = FALSE, ...),
      bw.f.BGM.boot1 = bw.f.BGM.boot1(y, w, kernel, ...),
      bw.f.BGM.boot2 = bw.f.BGM.boot2(y, w, kernel, plot = FALSE, ...),
      stop("unknown bandwidth rule")
    )
  }
  if (!is.finite(bw)) stop("non-finite 'bw'")
  if (bw <= 0) stop("'bw' is not positive")

  # Density Estimation
  aux <- bw^(-1) * outer(y.seq, y, "-")
  aux <- kernel_function_density(aux)
  aux <- aux %*% diag(weights)
  f.hat <- uw / bw * rowMeans(aux)

  if (plot == TRUE) {
    ord <- order(y.seq)
    plot(y.seq[ord], f.hat[ord],
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
    y.seq = y.seq,
    f.hat = f.hat,
    data.name = data.name,
    bw = bw,
    n = n,
    has.na = has.na,
    call = match.call()
  )
}
