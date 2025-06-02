#' Cross-validation bandwidth selection for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function performs bandwidth selection for \insertCite{jones1991;textual}{WData} kernel density estimation using cross-validation. It iterates through a range of bandwidth values and computes the cross-validation score for each bandwidth. The bandwidth that minimizes the cross-validation score is selected as the optimal bandwidth.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param lower Numeric value specifying the lower bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param upper Numeric value specifying the upper bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param nh An integer specifying the number of points for evaluating the cross-validation function. Default is 200.
#' @param tol Tolerance value used for checking if the minimum bandwidth occurs at the boundary of the range. Default is 10% of the lower bound.
#' @param plot Logical value indicating whether to plot the cross-validation function. Default is `TRUE`.
#' @return The optimal cross-validation bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator.
#' @details The optimal bandwidth is the one that minimizes the cross-validation function, i.e.,
#' \deqn{
#' \widehat{h}_{\mathrm{CV}}= \mathrm{arg} \min_{h>0} \mathrm{CV}(h) = \mathrm{arg} \min_{h>0} \int_{-\infty}^{+\infty} \widehat{f}_{\mathrm{J}}(y)^2 d y-2 \widehat{\mathbb{E}[\widehat{f}_{\mathrm{J}}]}.}
#' It holds that
#' \deqn{
#' \int_{-\infty}^{+\infty} \widehat{f}_{\mathrm{J}}(y)^2 dy = \frac{\widehat{\mu}_w^2}{n^{2} h} \sum_{i=1}^{n} \sum_{\substack{j=1}}^{n}
#' \frac{1}{w(Y_i)} \frac{1}{w(Y_j)}
#' (K \circ K) \left(\frac{Y_i-Y_j}{h} \right),
#' }
#' where \eqn{\circ} denotes convolution between two functions. Moreover,
#' \deqn{
#' \widehat{\mathbb{E}[\widehat{f}_{\mathrm{J}}]}=\frac{\widehat{\mu}_w}{n} \sum_{i=1}^n \frac{1}{w(Y_i)} \widehat{f}_{\mathrm{J}, -i}\left(Y_i\right),
#' \quad \text{where} \quad
#' \widehat{\mu}_w=n \left(\sum_{i=1}^{n}  \frac{1}{w(Y_i)} \right)^{-1}
#' }
#' and \eqn{\widehat{f}_{J, -i}} is \insertCite{jones1991;textual}{WData} estimator constructed without the \eqn{i}-th data point. In practice, it is not necessary to calculate \eqn{\widehat{f}_{J, -i}} for each data point in the sample, but the estimator
#' \eqn{\widehat{\mathbb{E}[\widehat{f}_{\mathrm{J}}]}} can be calculated as
#' \deqn{
#' \widehat{\mathbb{E}[\widehat{f}_{\mathrm{J}}]} = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n}  \frac{1}{w(Y_i)} \left( \sum_{j\neq i}  \frac{1}{w(Y_j)} \right)^{-1} \left(\sum_{j\neq i}  \frac{1}{w(Y_j)} K_h(Y_i-Y_j)\right).
#' }
#' The function optimizes the cross-validation function on the interval \eqn{I} determined by `lower` and `upper`. By default,
#' \eqn{I} is the one suggested by \insertCite{borrajo2017;textual}{WData}:
#' \deqn{
#' I = \left[\frac{\mathrm{IQR}}{2000n^{1/5}}, \frac{500 \mathrm{IQR} \log(n)^{1/5}}{n^{1/5}}\right],
#' }
#' where IQR is the interquartile range.
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @examples
#' bw.f.BGMcv(shrub.data$Width)
#' bw.f.BGMcv(shrub.data$Width, kernel = "epanechnikov")
#'
bw.f.BGMcv <-
  function(y,
           w = function(y) {
             ifelse(y >= 0, y, NA)
           },
           kernel = c(
             "gaussian", "epanechnikov",
             "rectangular", "triangular",
             "biweight", "cosine", "optcosine"
           ),
           lower = IQR(y) * n^(-0.2) * 2000^(-1),
           upper = IQR(y) * (log(n))^(0.2) * n^(-0.2) * 500,
           nh = 200L,
           tol = 0.1 * lower,
           plot = TRUE) {
    list2env(.check_biased_dataset(y, w), envir = environment())
    kernel <- match.arg(kernel)
    list2env(.get_kernel_values(kernel), envir = environment())

    hs <- .get_bandwidth_grid(nh, lower, upper, tol, plot)
    int <- numeric(nh)
    es <- numeric(nh)

    dist <- outer(y, y, FUN = "-")
    wi <- (sum(weights) - weights)^(-1)

    pb <- progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed, time to finish: :eta",
      total = nh,
      clear = FALSE,
      width = 60
    )

    for (k in 1:nh) {
      aux <- dist / hs[k]

      aux2 <- kernel_function_density(aux) %*% diag(weights) # (i,j) = K ((Yi - Yj)/h) * w(Yj) ^{-1}
      diag(aux2) <- 0
      es[k] <- -2 * uw * hs[k]^(-1) * mean(weights * wi * rowSums(aux2))

      aux3 <- diag(weights) %*% kernel_function_conv(aux) %*% diag(weights) # (i,j) = w(Yi) ^{-1} * w(Yj) ^{-1} * K o k((Yi - Yj)/h)
      int[k] <- uw^2 / n^2 / hs[k] * sum(aux3)

      pb$tick()
    }

    cvh <- int + es
    h <- hs[which.min(cvh)]
    if (h < lower + tol || h > upper - tol) {
      warning("minimum occurred at one end of the range")
    }

    if (plot == TRUE) {
      # Rule of the thumb with normal kernel as reference.

      sigma <- sqrt(uw * (mean(yw) - uw))
      bw.f.BGMnrd0 <- sigma * (8 * sqrt(pi) * RK * uw * uwb)^(0.2) * (3 * n * sigma_K_2^2)^(-0.2)

      # Plot the Cross - Validation function

      plot(
        hs, cvh,
        type = "l",
        xlab = "Bandwidth", ylab = "", sub = "CV(h)"
      )
      if (h < lower + tol || h > upper - tol) {
        title("Minimum occurred at one end of the range", col.main = "red")
      } else {
        title(paste0("CV bandwidth: ", round(h, 4)), col.main = "blue")
      }

      abline(v = h, h = min(cvh), col = "blue")
      abline(v = bw.f.BGMnrd0, lty = 2)

      par(ask = TRUE)

      plot(
        hs, int,
        type = "l",
        xlab = "Bandwidth",
        ylab = "",
        sub = "Bias(h)",
      )
      abline(v = h, col = "blue")
      abline(v = bw.f.BGMnrd0, lty = 2)

      plot(
        hs, es,
        type = "l",
        xlab = "Bandwidth",
        ylab = "",
        sub = "Var(h)",
      )
      abline(v = h, col = "blue")
      abline(v = bw.f.BGMnrd0, lty = 2)

      par(ask = FALSE)
    }

    return(h)
  }
