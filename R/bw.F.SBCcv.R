#' Cross-validation bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator
#'
#' This function performs bandwidth selection for \insertCite{bose2022;textual}{WData} kernel distribution estimator using cross-validation. It iterates through a range of bandwidth values and computes the cross-validation score for each bandwidth. The bandwidth that minimizes the cross-validation score is selected as the optimal bandwidth.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param lower Numeric value specifying the lower bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param upper Numeric value specifying the upper bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param nh An integer specifying the number of points for evaluation. Default is 200.
#' @param tol A numeric value specifying the tolerance for the bandwidth search. Default is \eqn{10\%} of the lower bound.
#' @param plot A logical value indicating whether to plot the cross-validation function. Default is `TRUE`.
#' @return The optimal bandwidth for \insertCite{bose2022;textual}{WData} kernel distribution estimator based on cross-validation.
#' @details The optimal bandwidth is obtained as the one that minimizes the cross-validation function, that is,
#' \deqn{\widehat{h}_{F, \mathrm{CV}} = \arg \min_{h_{F}>0} \frac{1}{n} \sum_{j=1}^n \left( \frac{\widehat{\mu}_w}{w(Y_j)} \mathbb{I} (y \geq Y_j) - \widehat{F}_{h_{F}, -j}(y)\right)^2 \!\!,
#' \quad \text{with} \quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}}
#' and \eqn{\widehat{F}_{h_{F}, -j}} is the \insertCite{bose2022;textual}{WData} kernel distribution estimator without the observation \eqn{Y_j}.
#' @references \insertAllCited{}
#' @seealso [`cdf.bd`][WData::cdf.bd()]
#' @examples
#' bw.F.SBCcv(shrub.data$Width)
#' bw.F.SBCcv(shrub.data$Width, kernel = "epanechnikov")
bw.F.SBCcv <- function(y,
                       w = function(y) {
                         ifelse(y >= 0, y, NA)
                       },
                       kernel = c(
                         "gaussian",
                         "epanechnikov",
                         "rectangular",
                         "triangular",
                         "biweight",
                         "cosine",
                         "optcosine"
                       ),
                       lower = IQR(y) * length(y)^(-1 / 3) * 0.05,
                       upper = IQR(y) * length(y)^(-1 / 3) * 5,
                       nh = 200L,
                       tol = 0.1 * lower,
                       plot = TRUE) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  hs <- .get_bandwidth_grid(nh, lower, upper, tol, plot)

  # Points for distribution evaluation
  from <- min(y) - (sort(y)[5] - min(y))
  to <- max(y) + (max(y) - sort(y, decreasing = T)[5])
  x <- seq.int(from, to, length.out = 511)
  uw_minus_j <- (n - 1) / (sum(weights) - weights)

  # Returns a matrix which stores by columns the estimations for each entry of "x" without one point.
  F_minus_j_h <- function(x, y, bw) {
    aux <- outer(x, y, "-") / bw # (i, j) = (x[i] - y[j]) / bw
    W_h_matrix <- kernel_function_distribution(aux) # (i, j) = N((x[i] - y[j]) / bw
    W_h_matrix_weights <- W_h_matrix %*% diag(weights) # (i,j) = N((x[i] - y[j]) / bw) * 1 / y[j]
    row_sums <- rowSums(W_h_matrix_weights)
    result <- ((row_sums - W_h_matrix_weights) %*% diag(uw_minus_j)) / (n - 1) # (i,j) = F_j,h(x[i])
    return(result)
  }

  cvh <- numeric(nh)

  # The estimation of the DF for each h.
  F_minus_j_h_vals <- lapply(hs, function(h) F_minus_j_h(x, y, h)) # A list of matrixs
  diff_matrix <- outer(x, y, "-") >= 0 # (i, j) = x[i] - y[j] >= 0
  weighted_diff <- diff_matrix %*% diag(weights) * uw # (i, j) = 1_{x[i] - y[j] >= 0} * 1 / y[j] *uw

  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed, time to finish: :eta",
    total = nh,
    clear = FALSE,
    width = 60
  )

  for (i in 1:nh) {
    integrand <- rowMeans((weighted_diff - F_minus_j_h_vals[[i]])^2)
    cvh[i] <- .simpsons_rule(x = x, fx = integrand)

    pb$tick()
  }

  h <- hs[which.min(cvh)]
  if (h < lower + tol || h > upper - tol) {
    warning("minimum occurred at one end of the range")
  }

  if (plot) {
    # Compute the rule of the thumb bandwidth
    sigma <- sqrt(uw * (mean(yw) - uw))
    bw.F.SBCnrd0 <- sigma * (sqrt(pi) * uw * uwb * kernel_kappa)^(1 / 3) * (n * kernel_eta^2)^(-1 / 3)

    plot(
      hs,
      cvh,
      type = "l",
      xlab = "Bandwidth",
      ylab = "",
      sub = "CV(h)"
    )
    if (h < lower + tol || h > upper - tol) {
      title("Minimum occurred at one end of the range", col.main = "red")
    } else {
      title(paste0("CV bandwidth: ", round(h, 4)), col.main = "blue")
    }

    abline(
      v = h,
      h = min(cvh),
      col = "blue"
    )
    abline(v = bw.F.SBCnrd0, lty = 2)
  }

  return(h)
}
