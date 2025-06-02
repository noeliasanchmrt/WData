#' Kernel quantile function estimator
#'
#' This function estimates the quantile function using the kernel method. It is a generalization of the quantile function estimator given in \insertCite{parzen1979;textual}{WData}  for unbiased data.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param x A numeric vector specifying the points where the quantile is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param bw The smoothing bandwidth to be used in the quantile estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, only option available is [`bw.q.SBCcv`][WData::bw.q.SBCcv()].
#' @param adjust A numeric value representing the manual adjustment factor for the bandwidth. Default is 1.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param from Numeric value specifying the lower bound for  evaluation when `x` is not provided. Default is 0.01.
#' @param to Numeric value specifying the upper bound for  evaluation when `x` is not provided. Default is 0.99.
#' @param nb An integer specifying the number of points for evaluation when `x` is not provided. Default is 99.
#' @param plot A logical value indicating whether to plot the quantile estimation. Default is `TRUE`.
#' @param ... Additional arguments to be passed to bandwidth selection functions.
#' @return A list with the following components:
#' \item{`x`}{The points where the quantiles are estimated.}
#' \item{`est_values`}{The estimated quantile values.}
#' \item{`bw`}{The bandwidth used.}
#' \item{`n`}{The sample size after elimination of missing values.}
#' \item{`call`}{The call which produced the result.}
#' \item{`data.name`}{The deparsed name of the y argument (biased dataset).}
#' \item{`has.na`}{Logical; indicates whether the original vector `y` contains any `NA` values.}
#' @details Estimator is given by
#' \deqn{\widehat{F_{h}^{-1}} (\tau) = \int_{0}^{1} \widehat{F^{-1}_n} (t) K_h \left(t-\tau\right) dt,}
#' where \eqn{\widehat{F^{-1}_n}} is the the empirical quantile function based on the estimator proposed by \insertCite{cox2005;textual}{WData} for the distribution function, \eqn{h} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h}(u)=1/h K\left(u / h\right)}.
#' @references \insertAllCited{}
#' @seealso [`qf.SBC`][WData::qf.SBC()]
#' @examples
#' qf.SBCkernel(shrub.data$Width, bw = 0.05, kernel = "epanechnikov")
#' plot(qf.SBC(y = shrub.data$Width), add = TRUE, col = "red") # For comparison
#'
qf.SBCkernel <- function(y,
                         w = function(y) {
                           ifelse(y >= 0, y, NA)
                         },
                         x,
                         bw = "bw.q.SBCcv",
                         adjust = 1,
                         kernel = c(
                           "gaussian", "epanechnikov",
                           "rectangular", "triangular",
                           "biweight", "cosine", "optcosine"
                         ),
                         from = 0.01,
                         to = 0.99,
                         nb = 99,
                         plot = TRUE,
                         ...) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  # Bandwidth
  if (is.character(bw)) {
    if (n < 2) stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(bw,
      bw.q.SBCcv = bw.q.SBCcv(y, w, kernel, plot = FALSE, ...),
      stop("unknown bandwidth rule")
    )
  }

  list2env(.get_prob_grid(x, from, to, nb, bw, adjust), envir = environment())

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw # Increasing sequence of quantiles
  # We force to have 10 points in each interval
  # t <- unlist(lapply(1:(length(ti) - 1), function(i) seq(ti[i], ti[i + 1], length.out = 26)))
  # t <- c(seq(0, ti[1], length.out = 26), t)
  # t <- sort(t)
  t <- seq(0, 1, length.out = 5000L)
  aux <- bw^(-1) * outer(x, t, "-")
  aux <- kernel_function_density(aux) / bw

  F_inv_hat <- matrix(0, nrow = length(x), ncol = length(t))
  for (j in 1:length(x)) {
    for (i in 1:n) {
      integral_value <- .simpsons_rule(
        t[(t >= ifelse(i == 1, 0, ti[i - 1])) & (t <= ifelse(i == n, 1, ti[i]))],
        aux[j, (t >= ifelse(i == 1, 0, ti[i - 1])) & (t <= ifelse(i == n, 1, ti[i]))]
      )
      F_inv_hat[j, i] <- vals[i] * integral_value
    }
  }

  yords <- rowSums(F_inv_hat) ## Estimations at x

  if (plot == TRUE) {
    order <- order(x)
    plot(x[order], yords[order],
      type = "l", main = "Quantile Function Estimator",
      xlab = paste(
        "n = ", n, "Bandwidth = ",
        format(round(bw, 5), nsmall = 5)
      ),
      ylab = "", col = "blue"
    )
    rug(ti)
  }

  list(
    x = x,
    est_values = yords,
    data.name = data.name,
    bw = bw,
    n = n,
    has.na = FALSE,
    call = match.call()
  )
}
