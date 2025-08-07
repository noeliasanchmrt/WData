#' \insertCite{efromovich2008;textual}{WData} adaptative density and distribution estimator
#'
#' This function computes \insertCite{efromovich2008;textual}{WData} density and distribution estimator given a biased dataset and its bias function.
#'
#' @param y  A numeric vector containing the biased data.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param x A numeric vector specifying the points where the distribution is estimated. Alternatively, `from`, `to` and `nb` can be used to define the evaluation points.
#' @param Jn Parameter which controls the maximum number of coefficients to be included in the low frequency part of the estimator.
#' @param Ct Parameters which controls the minimum size that the coefficients need to have to be included in the high frequency part of the estimation.
#' @param Cjm Parameter which controls the maximum number of coefficients to be included in the high frequency part of the estimator.
#' @param from Numeric value specifying the lower bound for  evaluation when `x` not provided. Default is calculated based on the range of input data.
#' @param to Numeric value specifying the upper bound for  evaluation when `x` not provided. Default is calculated based on the range of input data.
#' @param nb An integer specifying the number of points for evaluation when `x` not provided. Default is 512.
#' @param plot  Logical value indicating whether to plot the estimations. Default is `TRUE`.
#' @return A list with the following components:
#' \item{J}{The optimal cutoff \eqn{\widehat{J}}.}
#' \item{c}{The value for constant \eqn{c}.}
#' \item{U}{The value for the cutoff \eqn{U}.}
#' \item{x}{The points where the density is estimated.}
#' \item{fords}{The estimated density values.}
#' \item{Fords}{The estimated distribution values.}
#' \item{theta}{The estimated coefficients \eqn{\widehat{\theta_j}}.}
#' \item{omega}{The estimated weights \eqn{\widehat{\omega_j}}.}
#' \item{MISE}{MISE per each cutoff between zero and  `Cjm * Jn`.}
#' \item{x01}{The points in \eqn{[0,1]} where the density is estimated.}
#' \item{y01}{The estimator \eqn{\widehat{f}_{S, w, HF}(y)} evaluated in `x01`.}
#' @details Given a biased sample from a variable with support on \eqn{[0,1]},
#' the estimator proposed by \insertCite{efromovich2008;textual}{WData} is given by
#' \deqn{
#' \widehat{f}_{S, w, HF,c}(y) = \max \left( 0,  \widehat{f}_{S, w, HF}(y) - c\right),
#' }
#' where \eqn{c} is a constant ensuring that \eqn{\widehat{f}_{S, w, HF,c}}
#' integrates to \eqn{1} over \eqn{[0,1]}. Specifically,
#' \deqn{
#' \widehat{f}_{S, w, HF}(y) = \sum_{j=0}^{\widehat{J}} \widehat{\omega}_j \widehat{\theta}_j \varphi_j(y) + \sum_{j = \widehat{J}+1}^{c_{JM}J_n} \mathbb{I}\left(\widehat{\theta}_{j}^2 > c_T  \ln(n)/n\right) \widehat{\theta}_j \varphi _j (y)
#' }
#' where \eqn{\varphi_0(u)=1} and \eqn{\varphi_j(u)=\sqrt{2}\cos (\pi ju)} if
#' \eqn{j>0}. The first and second summations correspond to the low and high
#' frequency parts of the estimator, respectively. The parameters \eqn{c_{JM}}
#' and \eqn{J_n} control the maximum number of coefficients that constitute the
#' estimator \eqn{\widehat{f}_{S, w, HF}}. \insertCite{efromovich2008;textual}{WData}
#' suggests taking \eqn{c_{JM} = 6} and \eqn{c_{T}=4}. The weights
#' \eqn{\omega_j} are given by
#' \deqn{
#' \widehat{\omega}_{0} = 1
#' \quad \text{and} \quad
#' \widehat{\omega}_j = \max\left(0, 1 - (n \widehat{\theta}_j^2)^{-1}\right)
#' \quad j > 0
#' }
#' and the coefficients by
#' \deqn{
#' \widehat{\theta}_j =\widehat{\mu}_w n^{-1}\sum_{i=1}^{n} w^{-1}\left(Y_i\right) \varphi_j\left(Y_i\right) .
#' }
#' The threshold \eqn{\widehat{J}} is chosen as
#' \deqn{
#' \widehat{J} = \arg \min_{J \in [0,J_n]} \sum_{j=0}^{J}\left(2 n^{-1} - \widehat{\theta}_j^2 \right).
#' }
#' \insertCite{efromovich2008;textual}{WData} establishes \eqn{J_n = \lfloor 4 + 0.5
#' \ln{n} \rfloor} as a reasonable upper bound, where \eqn{\lfloor z \rfloor}
#' denotes the greatest integer less than or equal to \eqn{z}.
#'
#' The results presented for the interval \eqn{[0,1]} can be readily applied to
#' a compact interval \eqn{[a,b]}. For this purpose, given a sample \eqn{\{Y_1,
#' Y_2, \dots Y_n\}}, consider the transformation \eqn{Y_{i, [0,1]} = (Y_i -
#' a)/(b-a) \in [0,1]} and \eqn{w(x)_{[0,1]} = w((b-a)x+a)}, from which
#' \eqn{\widehat{f}_{S, w, HF,c, [0,1]}(x)} would be computed as indicated.
#' Finally, the estimation of the density function for the variable \eqn{X} is
#' given by
#' \deqn{
#' \widehat{f}_{S, w, HF, c, [a,b]}(y) = (b-a)^{-1} \widehat{f}_{S, w, HF, c, [0,1]}\left(\frac{y-a}{b-a}\right)
#' }
#' and its integral is taken as the estimation for the distribution function.
#' In many practical applications, the interval \eqn{[a,b]} is unknown. In such
#' cases, it is considered and interval of the form \eqn{[Y_{(1)} - \delta_1, Y_{(n)} +
#' \delta_2]}, where \eqn{Y_{(1)} \leq Y_{(2)} \leq \dots \leq Y_{(n)}} is the
#' ordered sample, and
#' \deqn{
#' \widehat{\delta_1} = Y_{(2)}- Y_{(1)}
#' \quad
#' \quad
#' \widehat{\delta_2} = Y_{(n)}- Y_{(n-1)}.
#' }
#' @references \insertAllCited{}
#' @examples
#' cdf.efro(y = shrub.data$Width)
cdf.efro <- function(y,
                     w = function(y) {
                       ifelse(y >= 0, y, NA)
                     },
                     x,
                     from,
                     to,
                     nb = 512L,
                     Jn = floor(4 + 0.5 * log(n)),
                     Ct = 4,
                     Cjm = 6,
                     plot = TRUE) {
  list2env(.check_biased_dataset(y, w), envir = environment())
  list2env(.get_xaxn_grid(y, x, from, to, nb, plot), envir = environment())

  message(sprintf(
    "Coefficient of difficulty due to biasing: %s",
    format(round(uw * uwb, 5), nsmall = 5)
  ))

  if (is.na(Jn) || Jn <= 0L || length(Jn) != 1) {
    stop("invalid Jn: Jn must be a positive number")
  }
  if (is.na(Ct) || Ct <= 0L || length(Ct) != 1) {
    stop("invalid Ct: Ct must be a positive number")
  }
  if (is.na(Cjm) || Cjm <= 0L || length(Cjm) != 1) {
    stop("invalid Cjm: Cjm must be a positive number")
  }

  y01 <- (y - from) / (to - from)

  phi <- function(u, j) {
    ifelse(j == 0, 1, sqrt(2) * cos(pi * j * u))
  }
  theta <- numeric(Cjm * Jn + 1)
  i <- 1
  for (i in 1:length(theta)) {
    theta[i] <- uw * n^{
      -1
    } * sum(weights * sapply(y01, FUN = phi, j = i))
    i <- i + 1
  }

  aux <- cumsum(2 / n - theta^2)
  J <- which.min(aux[1:(Jn + 1)])
  message(sprintf("Optimal Cutoff J: %d", J - 1))

  omega <- numeric(Cjm * Jn + 1)
  i <- 1
  for (i in 1:length(theta)) {
    omega[i] <- ifelse(i == 1, 1, max(0, 1 - theta[1]^2 / (n * theta[i]^
      2)))
    i <- i + 1
  }

  U <- Ct * theta[1] * log(n) / n
  message(sprintf(
    "Cutoff Value for High Frequency Terms U: %s",
    format(round(U, 5), nsmall = 5)
  ))
  HFT <- (theta^2 > U) * 1
  aux_wei <- c(omega[1:J], HFT[-c(1:J)])
  message(sprintf("Number of High Frequency Terms Included: %d", sum(HFT[-c(1:J)])))


  fSwHF01 <-
    function(x) {
      sum(aux_wei * theta * sapply(1:length(theta), FUN = phi, u = x))
    }


  x01 <- (x - from) / (to - from)
  yords01 <- sapply(x01, FUN = fSwHF01)

  c <- seq(min(yords01), max(yords01), len = 200)
  message(sprintf("Interval where c is searched: [%f, %f]", min(c), max(c)))
  cvalue <- numeric(length(c))

  i <- 1
  for (i in 1:length(c)) {
    f <- function(x) {
      max(c(0, fSwHF01(x = x) - c[i]), na.rm = TRUE)
    }
    cvalue[i] <- tryCatch(
      {
        integrate(Vectorize(f),
          lower = 0, upper = 1,
          subdivisions = 1000, rel.tol = .Machine$double.eps^.15
        )$value
      },
      error = function(e) {
        0
      }
    )

    i <- i + 1
  }

  copt <- c[which.min(abs(cvalue - 1))]
  cvalueopt <- cvalue[which.min(abs(cvalue - 1))]
  message(sprintf("Optimal Value for c: %s", format(round(copt, 5), nsmall = 5)))
  message(sprintf("Value of the Estimator Integral: %s", format(round(cvalueopt, 5), nsmall = 5)))


  yords01_positive <- pmax(0, yords01 - copt)

  yords <- (to - from)^
    {
      -1
    } * yords01_positive

  Fyords <- numeric(length(x))
  f <- approxfun(x, yords)
  for (i in 1:length(x)) {
    Fyords[i] <-
      integrate(f,
        lower = min(x), upper = x[i],
        rel.tol = .Machine$double.eps^0.21,
        subdivisions = 500L,
        stop.on.error = FALSE
      )$value
    i <- i + 1
  }

  if (plot == TRUE) {
    # par(mfrow = c(2, 4))
    plot(
      0:(length(theta) - 1),
      theta,
      type = "o", cex = 0.5,
      main = TeX(paste("$\\widehat{\\theta}$")),
      ylab = "", xlab = "Index",
      col = "grey"
    )
    abline(v = c(Jn, J - 1), lty = c(2, 3))
    abline(h = 0, col = "red")

    par(ask = TRUE)

    plot(
      0:(length(aux) - 1),
      aux,
      type = "o", cex = 0.5,
      main = "MISE",
      ylab = "", xlab = "Index",
      col = "grey",
      font.main = 1
    )
    abline(v = c(Jn, J - 1), lty = c(2, 3))
    plot(
      0:(length(omega) - 1),
      omega,
      type = "o", cex = 0.5,
      main = TeX(paste("$\\widehat{\\omega}$")),
      ylab = "", xlab = "Index",
      col = "grey"
    )
    abline(v = c(Jn, J - 1), lty = c(2, 3))
    plot(
      0:(length(theta) - 1),
      theta^2,
      type = "o", cex = 0.5,
      main = TeX(paste("$\\widehat{\\theta}^2$")),
      ylab = "", xlab = "Index",
      col = "grey"
    )
    abline(v = c(Jn, J - 1), lty = c(2, 3))
    abline(h = U, lty = 4)
    plot(
      x01,
      yords01,
      type = "l",
      main = TeX(paste("$\\widehat{f}_{S, w, HF, [0,1]}$")),
      ylab = "",
      xlab = "[0,1]"
    )
    abline(h = 0, col = "red")
    plot(
      x01,
      yords01_positive,
      type = "l",
      main = TeX(paste(
        "$\\max(0,\\widehat{f}_{S, w, HF, c, [0,1]} -c)$"
      )),
      ylab = "",
      xlab = "[0,1]"
    )
    abline(h = 0, col = "red")
    plot(
      x,
      yords,
      type = "l", col = "blue",
      main = TeX(paste("$\\widehat{f}_{S, w, HF, c, [a,b]}$")),
      ylab = "",
      xlab = paste0(
        "[",
        format(round(from, 2), nsmall = 2),
        ",",
        format(round(to, 2), nsmall = 2),
        "]"
      )
    )

    plot(
      x,
      Fyords,
      type = "l", col = "blue",
      main = TeX(paste("$\\widehat{F}_{S, w, HF, c, [a,b]}$")),
      ylab = "",
      xlab = paste0(
        "[",
        format(round(from, 2), nsmall = 2),
        ",",
        format(round(to, 2), nsmall = 2),
        "]"
      )
    )
    par(ask = FALSE)
  }

  return(invisible(list(
    J = J - 1,
    c = copt,
    U = U,
    x = x,
    fords = yords,
    Fords = Fyords,
    theta = theta,
    omega = omega,
    MISE = aux,
    x01 = x01,
    y01 = yords01
  )))
}
