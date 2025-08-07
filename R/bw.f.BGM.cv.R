#' Cross-validation bandwidth selector for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function estimates the bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator using cross-validation criteria from \insertCite{guillamon1998;textual}{WData}. It iterates through a range of bandwidth values and computes the cross-validation score for each bandwidth. The bandwidth that minimizes the cross-validation function is selected as the optimal bandwidth.
#'
#' @param y A numeric vector containing the biased sample.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample the sample `y`. By default, it is set to the length-biased function.
#' @param kernel A character string specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param lower Numeric value specifying the lower bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and sample size.
#' @param upper Numeric value specifying the upper bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and sample size.
#' @param nh An integer specifying the number of points to evaluate the cross-validation function. Default is 200.
#' @param tol Tolerance value used to check whether the minimum found lies at the boundaries of the interval; that is, the function will return a warning if the window minimizing the cross-validation function lies within `[lower, lower+tol]` or `[upper-tol, upper]`. Default is 10% of the lower bound.
#' @param plot Logical value indicating whether to plot the cross-validation function. Default is `TRUE`.
#' @return The optimal bandwidth value based on cross-validation criteria.
#' @details The optimal bandwidth is the one that minimizes the cross-validation function, i.e.,
#' \deqn{
#' \widehat{h}_{f, \mathrm{CV}}= \mathrm{arg} \min_{h_{f}>0} \mathrm{CV}(h_{f}) = \mathrm{arg} \min_{h_{f}>0} \int_{-\infty}^{+\infty} \widehat{f}_{\mathrm{J},h_{f}}(y)^2 d y-2 \widehat{\mathbb{E}}[\widehat{f}_{\mathrm{J},h_{f}}].}
#' It holds that
#' \deqn{
#' \int_{-\infty}^{+\infty} \widehat{f}_{\mathrm{J}, h_{f}}(y)^2 dy = \frac{\widehat{\mu}_w^2}{n^{2} h_{f}} \sum_{i=1}^{n} \sum_{\substack{j=1}}^{n}
#' \frac{1}{w(Y_i)} \frac{1}{w(Y_j)}
#' (K \circ K) \left(\frac{Y_i-Y_j}{h_{f}} \right),
#' }
#' where \eqn{\circ} denotes convolution between two functions and \eqn{\widehat{\mathbb{E}}[\widehat{f}_{\mathrm{J}, h_{f}}]} is computed as
#' \deqn{
#' \widehat{\mathbb{E}}[\widehat{f}_{\mathrm{J}, h_{f}}] = \frac{\widehat{\mu}_w}{n} \sum_{i=1}^{n}  \frac{1}{w(Y_i)} \left( \sum_{j\neq i}  \frac{1}{w(Y_j)} \right)^{-1} \left(\sum_{j\neq i}  \frac{1}{w(Y_j)} K_{h_{f}}(Y_i-Y_j)\right).
#' }
#' This function computes the bandwidth that minimizes the cross validation function, \eqn{\mathrm{CV}},  on the interval \eqn{I} determined by `lower` and `upper`. By default,
#' \eqn{I} is the one suggested by \insertCite{borrajo2017;textual}{WData}:
#' \deqn{
#' I = \left[\frac{\mathrm{IQR}}{2000n^{1/5}}, \frac{500 \mathrm{IQR} \log(n)^{1/5}}{n^{1/5}}\right],
#' }
#' where IQR is the interquartile range.
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @examples
#' bw.f.BGM.cv(shrub.data$Width)
#' bw.f.BGM.cv(shrub.data$Width, kernel = "epanechnikov")
#'
bw.f.BGM.cv <-
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
    list2env(.check_biased_sample(y, w), envir = environment())
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
      # Rule of thumb with gaussian kernel as reference.

      sigma <- sqrt(uw * (mean(yw) - uw))
      bw.f.BGM.rt <- sigma * (8 * sqrt(pi) * kernel_r * uw * uwb)^(0.2) * (3 * n * kernel_eta^2)^(-0.2)

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
      abline(v = bw.f.BGM.rt, lty = 2)

      par(ask = TRUE)

      plot(
        hs, int,
        type = "l",
        xlab = "Bandwidth",
        ylab = "",
        sub = "Bias(h)",
      )
      abline(v = h, col = "blue")
      abline(v = bw.f.BGM.rt, lty = 2)

      plot(
        hs, es,
        type = "l",
        xlab = "Bandwidth",
        ylab = "",
        sub = "Var(h)",
      )
      abline(v = h, col = "blue")
      abline(v = bw.f.BGM.rt, lty = 2)

      par(ask = FALSE)
    }

    return(h)
  }
