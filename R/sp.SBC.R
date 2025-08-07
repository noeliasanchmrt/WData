#' Kernel sparsity estimator
#'
#' This function estimates the sparsity function using the kernel method. It is a generalization of the quantile function estimator given in \insertCite{parzen1979;textual}{WData} and \insertCite{falk1986;textual}{WData} for unbiased data.
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
#' @details The sparsity estimator is given by:
#' \deqn{\widehat{\left(F^{-1}\right)_{h}^{(1)}} (\tau) = \sum_{i=1}^{n} Y_{(i)} \left( K_h \left( T_{i-1} - \tau \right) - K_h \left( T_{i} - \tau \right) \right),
#' \quad \text{where} \quad
#' T_i = \sum_{j=1}^{i} \frac{1}{Y_{(j)}} \bigg/ \sum_{i=j}^{n} \frac{1}{Y_{(j)}}.}
#' where \eqn{Y_{(i)}} is the \eqn{i}-th order statistic of the sample, \eqn{h} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h}(u)=1/h K\left(u / h\right)}.
#' @references \insertAllCited{}
#' @seealso [`qf.SBC`][WData::qf.SBC()]
#' @examples
#' sp.SBC(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#'
sp.SBC <- function(y,
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
  aux <- bw^(-1) * outer(x, c(0, ti), "-") # We include the T_0 = 0
  aux <- kernel_function_density(aux)
  diff <- c(vals[1], diff(vals), -vals[n]) # We include Y(1) and -Y(n)
  product <- aux %*% diag(diff)
  yords <- 1 / bw * rowSums(product)

  if (plot == TRUE) {
    order <- order(x)
    plot(x[order], yords[order],
      type = "l", main = "SBC's Sparsity Estimator",
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

    diff <- c(vals[1], diff(vals), vals[n]) # We include Y(1) and -Y(n)
    product <- aux %*% diag(diff)
    yords <- 1 / bw * rowSums(product)
    # lines(x, yords, col = "red")

    diff <- c(vals[1], diff(vals), -vals[n]) # We include Y(1) and -Y(n)
    product <- aux %*% diag(diff)
    product <- aux[, 2:(n - 1)] %*% diag(diff[2:(n - 1)])
    yords <- 1 / bw * rowSums(product)
    # lines(x, yords, col = "magenta")
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
