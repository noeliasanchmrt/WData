#' Cross-validation bandwidth selection for kernel quantile estimator
#'
#' This function performs bandwidth selection for kernel quantile estimation using cross-validation. It iterates through a range of bandwidth values and computes the cross-validation score for each bandwidth. The bandwidth that minimizes the cross-validation score is selected as the optimal bandwidth.
#'
#' @useDynLib WData
#' @importFrom Rcpp evalCpp
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param lower Numeric value specifying the lower bound for bandwidth selection. Default is 0.0001.
#' @param upper Numeric value specifying the upper bound for bandwidth selection. Default is 0.2.
#' @param nh An integer specifying the number of points for evaluating the cross-validation function. Default is 200.
#' @param tol Tolerance value used for checking if the minimum bandwidth occurs at the boundary of the range. Default is 10% of the lower bound.
#' @param plot Logical value indicating whether to plot the cross-validation function. Default is `TRUE`.
#' @return The optimal cross-validation bandwidth for kernel quantile estimator.
#' @details The optimal bandwidth is the one that minimizes the cross-validation function. Specifically, the bandwidth \eqn{\widehat{h}_{\operatorname{CV}}} is defined as
#' \deqn{
#'   \widehat{h}_{\operatorname{CV}} = \arg \min_{h \in [h_1, h_2]} \frac{1}{n} \sum_{i=1}^{n} \left(Y_{(i)} - \widehat{F_{h,-i}^{-1}}(T_i)\right)^2,
#' }
#' where \eqn{\widehat{F_{h,-i}^{-1}}} denotes the quantile function estimator computed without the \eqn{i}-th order statistic.
#' @references \insertAllCited{}
#' @seealso [`qf.SBCkernel`][WData::qf.SBCkernel()]
#' @examples
#' bw.q.SBCcv(shrub.data$Width)
#' bw.q.SBCcv(shrub.data$Width, kernel = "epanechnikov")
#'
bw.q.SBCcv <- function(y, w = function(y) {
                         ifelse(y >= 0, y, NA)
                       }, kernel = c(
                         "gaussian", "epanechnikov", "rectangular", "triangular",
                         "biweight", "cosine", "optcosine"
                       ),
                       lower = 0.0001,
                       upper = 0.2,
                       nh = 200L,
                       tol = 0.1 * lower,
                       plot = TRUE) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  hs <- .get_bandwidth_grid(nh, lower, upper, tol, plot)
  cvh <- numeric(nh)

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw # Ti when sample is complete
  Tis <- t(sapply(1:n, function(i) {
    w <- weightsvals[-i]
    cumsum(w) / sum(w)
  })) # (i,j) = T(j) when (i) is not in the sample

  vals_minus_i <- t(sapply(1:n, function(i) {
    vals[-i]
  })) # (i,j) = Y(j) when (i) is not in the sample


  start <- Sys.time()

  cvh <- .Call("_WData_quantile_cross_validation", Tis, vals_minus_i, ti, vals, hs, kernel)

  end <- Sys.time()
  message(sprintf(
    "Execution time in C++: %s",
    format(round(end - start, 5), nsmall = 5)
  ))

  h <- hs[which.min(cvh)]
  if (h < lower + tol || h > upper - tol) {
    warning("minimum occurred at one end of the range")
  }
  if (plot == TRUE) {
    plot(hs, cvh,
      type = "l", xlab = "Bandwidth", ylab = "",
      sub = "CV(h)"
    )
    if (h < lower + tol || h > upper - tol) {
      title("Minimum occurred at one end of the range",
        col.main = "red"
      )
    } else {
      title(paste0("Optimal Bandwidth h:", round(h, 4)),
        col.main = "blue"
      )
    }
    abline(v = h, h = min(cvh), col = "blue")
  }
  return(h)
}
