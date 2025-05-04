#' Quantile-quantile plots
#'
#' This function produces a normal Q-Q plot of the values in `y` with the bias function `w`.
#'
#' @param y The biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param ylim The limits for the axis.
#' @param main An overall title for the plot.
#' @param xlab A label for the x axis.
#' @param ylab A label for the y axis.
#' @param plot.it A logical value indicating whether to plot the Q-Q plot. Default is `TRUE`.
#' @param datax A logical value indicating whether the data should be plotted on the x-axis. Default is `FALSE`.
#' @param probs Numeric vector of length two, representing probabilities. Corresponding quantile pairs define the line drawn.
#' @param distribution Quantile function for reference theoretical distribution.
#' @param qmethod Character string specifying the quantile estimator to be used. Options are `"qf.SBC"` ([`qf.SBC`][WData::qf.SBC()]) and `"qf.sen"` ([`qf.sen`][WData::qf.sen()]). Default is `"qf.SBC"`.
#' @param ... Additional graphical parameters.
#' @return For `w_qqnorm`, a list with components:
#' \item{x}{The `x` coordinates of the points that were/would be plotted}
#' \item{y}{The original `y` vector, i.e., the corresponding `y` coordinates including `NA`s.}
#' @details
#' This function is a modification of the [`qqnorm`][stats::qqnorm()] function, which is used to create boxplots for biased data. [`w_qqnorm`][WData::w_qqnorm()] is a generic function whose default method produces a normal Q-Q plot of the values in `y`. [`w_qqline`][WData::w_qqline()] adds a line to a theoretical Q-Q plot, by default normal, which passes through the quantiles defined in `probs`, by default the first and third quartiles.
#' @seealso [`qf.SBC`][WData::qf.SBC()], [`qf.sen`][WData::qf.sen()]
#' @examples
#' w_qqnorm(shrub.data$Width)
#' w_qqline(shrub.data$Width, col = "blue", qmethod = "qf.sen")
#' w_qqline(shrub.data$Width, col = "magenta", qmethod = "qf.SBC")
#' @aliases w_qqnorm w_qqnorm.default w_qqline
#' @name w_qqnorm

#' @rdname w_qqnorm
#' @order 1
#' @export w_qqnorm
w_qqnorm <- function(y, ...) {
  UseMethod("w_qqnorm")
}

#' @rdname w_qqnorm
#' @order 2
#' @export w_qqnorm.default
w_qqnorm.default <-
  function(y, w = function(y) {
             ifelse(y >= 0, y, NA)
           }, ylim, main = "Normal Q-Q Plot",
           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
           plot.it = TRUE, datax = FALSE, ...) {
    if (has.na <- any(ina <- is.na(y))) { ## keep NA's in proper places
      yN <- y
      y <- y[!ina]
    }
    if (0 == (n <- length(y))) {
      stop("y is empty or has only NAs")
    }

    list2env(.check_biased_dataset(y, w), envir = environment())

    if (plot.it && missing(ylim)) ylim <- range(y)

    vals <- sort(y)
    weightsvals <- sapply(vals, w)^(-1)
    ti <- cumsum(weightsvals) / n * uw
    ti <- (c(0, ti) + c(ti, 1)) / 2 # Adjustment! In order to match ppoints(n)

    ti <- ifelse(ti >= 1, 0.99, ti) # Adjustment for the qnorm function
    x <- qnorm(ti)[order(order(y))]

    if (has.na) {
      y <- x
      x <- yN
      x[!ina] <- y
      y <- yN
    }

    if (plot.it) {
      if (datax) {
        plot(y, x, main = main, xlab = ylab, ylab = xlab, xlim = ylim, ...)
      } else {
        plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
      }
    }
    invisible(if (datax) list(x = y, y = x) else list(x = x, y = y))
  }

#' @rdname w_qqnorm
#' @order 3
#' @export w_qqline
w_qqline <- function(y, w = function(y) {
                       ifelse(y >= 0, y, NA)
                     }, datax = FALSE, distribution = qnorm,
                     probs = c(0.25, 0.75), qmethod = c("qf.SBC", "qf.sen"), ...) {
  stopifnot(length(probs) == 2, is.function(distribution))

  # Since we are calling either qf.SBC or qf.sen we don't need to check y and w again.
  qmethod <- match.arg(qmethod)

  y <- switch(qmethod,
    qf.SBC = qf.SBC(y, w = w)(probs),
    qf.sen = qf.sen(y, w = w)(probs),
    stop("unknown method for quantile computation")
  )

  x <- distribution(probs)
  if (datax) {
    slope <- diff(x) / diff(y)
    int <- x[1L] - slope * y[1L]
  } else {
    slope <- diff(y) / diff(x)
    int <- y[1L] - slope * x[1L]
  }
  abline(int, slope, ...)
}
