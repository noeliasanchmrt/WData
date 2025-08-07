#' Create an "eqf" structure
#'
#' This function generates an object of class `"eqf"` containing a vector and its associated probabilities.
#'
#' @param x A numeric vector or other data type that serves as the basis for the structure.
#' @param probs A numeric vector of probabilities (default is `seq(0, 1, 0.25)`).
#' @param xlab,ylab A character string for the x-axis/y-axis label for method [`plot`][graphics::plot()].
#' @param verticals A logical value indicating whether to add vertical lines to the plot (default is `FALSE`).
#' @param col.01line A character string for the color of the 0-1 lines (default is `"gray70"`).
#' @param pch An integer specifying the plotting character for method [`plot`][graphics::plot()] (default is 19).
#' @param digits An integer specifying the number of significant digits to be printed (default is `getOption("digits") - 2L`).
#' @param object An object of class `"eqf"`.
#' @param ... Additional arguments to be passed to methods.
#' @return An object of class `"eqf"`, which is a list with elements `x` and `probs`.
#' @seealso [`qf.sen`][WData::qf.sen()], [`qf.SBC`][WData::qf.SBC()]
#' @examples
#' eqf(1:5, (1:5) / 5)
#' @aliases eqf plot.eqf print.eqf summary.eqf quantile.eqf
#' @name eqf
#' @export

#' @rdname eqf
#' @order 1
#' @export eqf
eqf <- function(x, probs = seq(0, 1, 0.25)) {
  x <- sort(unique(x))
  probs <- sort(probs)

  eqf_fun <- stepfun(probs, c(0, x))

  attr(eqf_fun, "x") <- x
  attr(eqf_fun, "probs") <- probs
  class(eqf_fun) <- c("eqf", class(eqf_fun))

  return(eqf_fun)
}

#' @rdname eqf
#' @order 2
#' @export plot.eqf
plot.eqf <- function(x, ..., ylab = "", xlab = "", verticals = FALSE, col.01line = "gray70",
                     pch = 19) {
  plot.stepfun(x, ...,
    ylab = ylab, verticals = verticals,
    pch = pch
  )
  abline(v = c(0, 1), col = col.01line, lty = 2)
}

#' @rdname eqf
#' @order 3
#' @export print.eqf
print.eqf <- function(x, digits = getOption("digits") - 2L, ...) {
  numform <- function(x) {
    paste(formatC(x, digits = digits),
      collapse = ", "
    )
  }
  cat("Empirical QF \nCall: ")
  print(attr(x, "call"), ...)
  n <- length(xx <- environment(x)$x)
  i1 <- 1L:min(3L, n)
  i2 <- if (n >= 4L) {
    max(4L, n - 1L):n
  } else {
    integer()
  }
  cat(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3L) {
    ", "
  }, if (n > 5L) {
    " ..., "
  }, numform(xx[i2]), "\n", sep = "")
  invisible(x)
}

#' @rdname eqf
#' @order 4
#' @export summary.eqf
summary.eqf <- function(object, ...) {
  n <- length(eval(expression(x), envir = environment(object)))
  header <- paste("Empirical CDF:\t ", n, "unique values with summary\n")
  structure(summary(knots(object), ...), header = header, class = "summary.ecdf")
}

#' @rdname eqf
#' @order 5
#' @export quantile.eqf
quantile.eqf <- function(x, ...) {
  quantile(
    evalq(rep.int(x, diff(c(0, round(nobs * y)))), environment(x)),
    ...
  )
}
