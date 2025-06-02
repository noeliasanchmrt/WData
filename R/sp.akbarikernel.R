#' \insertCite{akbari2019;textual}{WData} kernel sparsity estimator
#'
#' This function estimates the sparsity of a distribution using \insertCite{akbari2019;textual}{WData} kernel estimator.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param x A numeric vector specifying the points where the sparsity is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The smoothing bandwidth to be used in the sparsity estimation.
#' @param adjust A numeric value representing the manual adjustment factor for the bandwidth. Default is 1.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound for  evaluation when `x` is not provided. Default is 0.01.
#' @param to Numeric value specifying the upper bound for  evaluation when `x` is not provided. Default is 0.99.
#' @param nb An integer specifying the number of points for evaluation when `x` is not provided. Default is 99.
#' @param plot A logical value indicating whether to plot the sparsity estimation. Default is `TRUE`.
#' @return A list with the following components:
#' \item{`x`}{The points where the sparsity is estimated.}
#' \item{`est_values`}{The estimated sparsity values.}
#' \item{`bw`}{The bandwidth used.}
#' \item{`n`}{The sample size after elimination of missing values.}
#' \item{`call`}{The call which produced the result.}
#' \item{`data.name`}{The deparsed name of the y argument (biased dataset).}
#' \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NA` values.}
#' @details \insertCite{akbari2019;textual}{WData} sparsity estimator is given by
#' \deqn{\widehat{\left(F^{-1}\right)_{h}^{(1)}} (\tau) = \int_{0}^{1} K_h (t-\tau) d \widehat{F^{-1}_{\operatorname{SEN}}} (t),}
#' where \eqn{\widehat{F^{-1}_{\operatorname{SEN}}}} is \insertCite{sen1984;textual}{WData} quantile estimator, \eqn{h} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h}(u)=1/h K\left(u / h\right)}.
#' @references \insertAllCited{}
#' @seealso [`qf.sen`][WData::qf.sen()]
#' @examples
#' sp.akbarikernel(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#'
sp.akbarikernel <- function(y,
                            w = function(y) {
                              ifelse(y >= 0, y, NA)
                            },
                            x,
                            bw = NULL,
                            adjust = 1,
                            kernel = c(
                              "gaussian", "epanechnikov",
                              "rectangular", "triangular",
                              "biweight", "cosine", "optcosine"
                            ),
                            from = 0.01,
                            to = 0.99,
                            nb = 99,
                            plot = TRUE) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  list2env(.get_prob_grid(x, from, to, nb, bw, adjust), envir = environment())

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw # Increasing sequence of quantiles
  aux <- bw^(-1) * outer(x, ti, "-")
  aux <- kernel_function_density(aux)
  diff <- c(vals[1], diff(vals)) # (Y_(i)- Y(i-1))
  aux <- aux %*% diag(diff)
  yords <- 1 / bw * rowSums(aux)

  if (plot == TRUE) {
    order <- order(x)
    plot(x[order], yords[order],
      type = "l", main = "Akbari's Sparsity Estimator",
      xlab = paste(
        "n = ", n, "Bandwidth = ",
        format(round(bw, 5), nsmall = 5)
      ),
      ylab = "Density", col = "blue"
    )
    rug(ti)
    legend("topright",
      legend = c("Estimated Sparsity"),
      col = c("blue"),
      lty = c(1), cex = 0.8
    )
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
