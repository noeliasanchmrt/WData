#' \insertCite{akbari2019;textual}{WData} plug-in sparsity estimator
#'
#' This function estimates the sparsity of a distribution using \insertCite{akbari2019;textual}{WData} plug-in sparsity estimator.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param bw The smoothing bandwidth to be used in the density estimation. `bw` can also be a character string giving a rule to choose the bandwidth. In this case, options available are [`bw.f.BGMnrd0`][WData::bw.f.BGMnrd0()], [`bw.f.BGMcv`][WData::bw.f.BGMcv()],  [`bw.f.BGMboot1`][WData::bw.f.BGMboot1()] and [`bw.f.BGMboot2`][WData::bw.f.BGMboot2()].
#' Default is [`bw.f.BGMnrd0`][WData::bw.f.BGMnrd0()]. The specified (or computed) value of `bw` is multiplied by `adjust`.
#' @param adjust A numeric value representing the manual adjustment factor for the bandwidth. Default is 1.
#' @param kernel A character vector specifying the kernel function. Available options: `"gaussian"`, `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` and `"optcosine"`.
#' @return A function of class `eqf`, inheriting from the [`stepfun`][stats::stepfun()] class, and hence inheriting a [`knots`][stats::knots()] method.
#' @details \insertCite{akbari2019;textual}{WData} suggests estimating the sparsity function using the following plug-in estimator:
#' \deqn{\widehat{\left(F^{-1}\right)_{\operatorname{PI}, h}^{(1)}} (\tau) =
#' \frac{1}{\widehat{f}_{\operatorname{J}, h}\left( \widehat{F^{-1}_{\operatorname{SEN}}} (\tau )\right)},}
#' where \eqn{\widehat{f}_{\operatorname{J}, h}} is \insertCite{jones1991;textual}{WData} kernel density estimator and \eqn{\widehat{F^{-1}_{\operatorname{SEN}}}} is \insertCite{sen1984;textual}{WData} quantile function estimator.
#' @references \insertAllCited{}
#' @seealso [`qf.sen`][WData::qf.sen()], [`df.jones`][WData::df.jones()]
#' @examples
#' sp.akbaripi(shrub.data$Width, bw = 0.25)
sp.akbaripi <- function(y,
                        w = function(y) {
                          ifelse(y >= 0, y, NA)
                        },
                        bw = "bw.f.BGMnrd0",
                        adjust = 1,
                        kernel = c(
                          "gaussian", "epanechnikov",
                          "rectangular", "triangular",
                          "biweight", "cosine", "optcosine"
                        )) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  kernel <- match.arg(kernel)
  list2env(.get_kernel_values(kernel), envir = environment())

  vals <- sort(y)
  weightsvals <- sapply(vals, w)^(-1)
  ti <- cumsum(weightsvals) / n * uw
  jones <- df.jones(y = vals, w = w, bw = bw, adjust = adjust, kernel = kernel, x = vals, plot = FALSE)
  rval <- approxfun(ti, 1 / jones$est_values,
    method = "constant", rule = 2,
    f = 0, ties = "ordered"
  )
  class(rval) <- c("eqf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
