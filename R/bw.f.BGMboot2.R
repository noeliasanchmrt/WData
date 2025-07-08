#' \insertCite{borrajo2017;textual}{WData} bootstrap bandwidth selection for \insertCite{jones1991;textual}{WData} kernel density estimator
#'
#' This function estimates the bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator using the bias-corrected bootstrap method developed by \insertCite{borrajo2017;textual}{WData} and based on the methodology introduced by \insertCite{bose2013;textual}{WData}.
#'
#' @param y A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @param bw0 A character string specifying the pilot bandwidth. It can also be a character string giving a rule to choose the bandwidth. Options available can be checked in [`bw.nrd`][stats::bw.nrd()]. By default it is set to `1/8 n^(-1/9)`.
#' @param lower Numeric value specifying the lower bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param upper Numeric value specifying the upper bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and number of data points.
#' @param nh An integer specifying the number of points for evaluating the MISE function. Default is 200.
#' @param tol  Tolerance value used for checking if the minimum bandwidth occurs at the boundary of the range. Default is 10% of the lower bound.
#' @param from Numeric value specifying the lower bound for  evaluation. Default is calculated based on the range of input data.
#' @param to Numeric value specifying the upper bound for  evaluation. Default is calculated based on the range of input data.
#' @param plot  Logical value indicating whether to plot the function to be minimized. Default is `TRUE`.
#' @return The estimated bootstrap bandwidth for \insertCite{jones1991;textual}{WData} kernel density estimator.
#' @details
#' The bandwidth returned is the one minimizing \eqn{\mathrm{MISE}^{\ast}} over a compact interval \eqn{[h_1,h_2]} (determined by arguments `lower` and `upper`), i.e.,
#' \deqn{
#' \widehat{h}_{f, \mathrm{B}} = \arg \min_{h_{f} \in [h_1,h_2]} \int_{-\infty}^{+ \infty} \mathrm{MSE}^{\ast}\left(\widehat{f}^{\ast}_{\mathrm{J}, h_{f}}(y)\right)dy.
#' }
#' \eqn{\mathrm{MISE}^{\ast}} and \eqn{\mathrm{MSE}^{\ast}} correspond with the expression of the mean integrated squared error and the mean squared error of the bootstrap estimator \eqn{\widehat{f}^{\ast}_{\mathrm{J}, h_{f}}} provided by \insertCite{borrajo2017;textual}{WData}.
#' @references \insertAllCited{}
#' @seealso [`df.jones`][WData::df.jones()]
#' @examples
#' bw.f.BGMboot2(shrub.data$Width, nh = 50L)
bw.f.BGMboot2 <- function(y,
                          w = function(y) {
                            ifelse(y >= 0, y, NA)
                          },
                          kernel = c(
                            "gaussian", "epanechnikov",
                            "rectangular", "triangular",
                            "biweight", "cosine", "optcosine"
                          ),
                          bw0 = 1 / 8 * n^(-1 / 9),
                          lower = IQR(y) * n^(-0.2) * 2000^(-1),
                          upper = IQR(y) * (log(n) / n)^(0.2) * 500,
                          nh = 200L,
                          tol = 0.1 * lower,
                          from = min(y) - (sort(y)[5] - min(y)),
                          to = max(y) + (max(y) - sort(y, decreasing = T)[5]),
                          plot = TRUE) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())
  hs <- .get_bandwidth_grid(nh, lower, upper, tol, plot)

  if (!is.numeric(from) || is.na(from) || length(from) != 1) {
    stop("invalid from: from must be numeric")
  }
  if (!is.numeric(to) || is.na(to) || length(to) != 1) {
    stop("invalid to: to must be numeric")
  }
  if (from >= to) {
    stop("from must be less than to")
  }

  message(sprintf("Interval where density is evaluated: [%f, %f]", from, to))

  if (!is.finite(bw0)) stop("non-finite 'bw0'")
  if (bw0 <= 0) stop("'bw0' is not positive")

  # We calculate a biased density estimation (hat_gnh0)
  gnbw0 <- density(
    x = y, from = from, to = to,
    n = 511, bw = bw0, kernel = kernel
  )
  message(sprintf("Pilot Bandwidth for Bootstrap: %f", gnbw0$bw))

  wz <- 1 / sapply(gnbw0$x, w)
  if (any(is.na(wz))) {
    warning("w is not evaluable at some points of the interval [from, to]")
  }
  index1 <- !is.na(wz)

  wz <- wz[index1]
  gnbw0$x <- gnbw0$x[index1]
  gnbw0$y <- gnbw0$y[index1]

  gnbw0wz <- gnbw0$y * wz

  # aux1 <- cubintegrate(approxfun(gnbw0$x, gnbw0wz), lower = from, upper = to)$integral
  aux1 <- .simpsons_rule(gnbw0$x, gnbw0wz)

  Bias <- Variance <- numeric(length(hs))

  fjh <- list()
  fjh$x <- seq.int(from = from, to = to, length.out = 511)
  fjh$y <- numeric(length(fjh$x))

  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed, time to finish: :eta",
    total = nh,
    clear = FALSE,
    width = 60
  )

  for (k in 1:nh) {
    aux <- outer(fjh$x, y, "-") / hs[k]
    aux <- kernel_function_density(aux)
    aux <- aux %*% diag(weights)
    fjh$y <- uw * rowMeans(aux) / hs[k]

    if (any(is.na(fjh$y))) {
      warning("'Jones' estimation could not be computed for some bandwidth")
      index2 <- !is.na(fjh$y)
      fjh$x <- fjh$x[index2]
      fjh$y <- fjh$y[index2]
    }

    aux2 <- aux3 <- numeric(length(fjh$x))

    aux_khyz <- outer(fjh$x, gnbw0$x, "-") / hs[k]
    Khyz <- kernel_function_density(aux_khyz) / hs[k]
    gnbw0wzKhyz <- (Khyz %*% diag(wz)) %*% diag(gnbw0$y)
    gnbw0wzKhyz2 <- (Khyz %*% diag(wz))^2 %*% diag(gnbw0$y)

    # for (j in 1:length(fjh$x)) {
    #   aux2[j] <- cubintegrate(approxfun(gnbw0$x, gnbw0wzKhyz[j, ], rule = 2),
    #     lower = from, upper = to
    #   )$integral
    #   aux3[j] <- cubintegrate(approxfun(gnbw0$x, gnbw0wzKhyz2[j, ], rule = 2),
    #     lower = from, upper = to
    #   )$integral
    # }

    for (j in 1:length(fjh$x)) {
      aux2[j] <- .simpsons_rule(gnbw0$x, gnbw0wzKhyz[j, ])
      aux3[j] <- .simpsons_rule(gnbw0$x, gnbw0wzKhyz2[j, ])
    }

    aux2aux1dfhy2 <- (aux2 / aux1 - fjh$y)^2
    if (any(is.na(aux2aux1dfhy2))) {
      warning("potencial error in the calculations")
    }
    aux3aux22 <- n^(-1) * aux1^(-2) * (aux3 - (2 * n - 1) * aux2^2) # Corrected (not change)
    if (any(is.na(aux3aux22))) {
      warning("potencial error in the calculations")
    }

    # Bias[k] <- cubintegrate(approxfun(fjh$x, aux2aux1dfhy2, rule = 2),
    #   lower = from, upper = to
    # )$integral
    # Variance[k] <- cubintegrate(approxfun(fjh$x, aux3aux22, rule = 2),
    #   lower = from, upper = to
    # )$integral

    Bias[k] <- .simpsons_rule(fjh$x, aux2aux1dfhy2)
    Variance[k] <- .simpsons_rule(fjh$x, aux3aux22)

    pb$tick()
  }


  if (any(is.na(Bias))) {
    warning("bias could not be computed for some bandwidth")
  }

  if (any(is.na(Variance))) {
    warning("variance could not be computed for some bandwidth")
  }

  MISE <- Bias + Variance
  paste(c("MISE: ", MISE), collapse = " ")

  if (any(is.na(MISE))) {
    warning("MISE could not be computed for some bandwidth")
  }

  if (all(is.na(MISE))) {
    stop("MISE could not be computed for any bandwidth")
  }

  h <- hs[which.min(MISE)]

  if (!identical(h, numeric(0))) {
    if (h < lower + tol || h > upper - tol) {
      warning("minimum occurred at one end of the range")
    }
  }

  h <- ifelse(!identical(h, numeric(0)), h, NA)

  if (plot == TRUE) {
    # par(mfrow = c(1, 3))
    plot(
      hs, MISE,
      type = "l", xlab = "Bandwidth", ylab = "", sub = "MISE(h)"
    )

    abline(v = h, h = min(MISE), lty = 1, col = "blue")
    abline(v = bw.f.BGMnrd0(y, w, kernel = kernel), lty = 2)
    rug(hs)

    if (h < lower + tol || h > upper - tol) {
      title("Minimum occurred at one end of the range", col.main = "red")
    } else {
      title(paste0("Bootstrap bandwidth: ", round(h, 4)), col.main = "blue")
    }

    plot(hs, Bias,
      type = "l", xlab = "Bandwidth", ylab = "", sub = "Bias(h)"
    )
    abline(v = h, h = min(Bias), lty = 1, col = "blue")
    abline(v = bw.f.BGMnrd0(y, w, kernel = kernel), lty = 2)
    rug(hs)

    par(ask = TRUE)

    plot(hs, Variance,
      type = "l", xlab = "Bandwidth", ylab = "", sub = "Var(h)"
    )
    abline(v = h, h = min(Variance), lty = 1, col = "blue")
    abline(v = bw.f.BGMnrd0(y, w, kernel = kernel), lty = 2)
    rug(hs)

    par(ask = FALSE)
  }

  return(h)
}
