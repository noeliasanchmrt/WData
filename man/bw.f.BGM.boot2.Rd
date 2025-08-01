% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bw.f.BGM.boot2.R
\name{bw.f.BGM.boot2}
\alias{bw.f.BGM.boot2}
\title{\insertCite{borrajo2017;textual}{WData} bootstrap bandwidth selector for \insertCite{jones1991;textual}{WData} kernel density estimator}
\usage{
bw.f.BGM.boot2(
  y,
  w = function(y) {
     ifelse(y >= 0, y, NA)
 },
  kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight",
    "cosine", "optcosine"),
  bw0 = 1/8 * n^(-1/9),
  lower = IQR(y) * n^(-0.2) * 2000^(-1),
  upper = IQR(y) * (log(n)/n)^(0.2) * 500,
  nh = 200L,
  tol = 0.1 * lower,
  from = min(y) - (sort(y)[5] - min(y)),
  to = max(y) + (max(y) - sort(y, decreasing = TRUE)[5]),
  plot = TRUE
)
}
\arguments{
\item{y}{A numeric vector containing the biased sample.}

\item{w}{A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample \code{y}. By default, it is set to the length-biased function.}

\item{kernel}{A character string specifying the kernel function. Available options: \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, \code{"biweight"}, \code{"cosine"} and \code{"optcosine"}.}

\item{bw0}{The bandwidth value to be used in \code{\link[stats:density]{density}}. It can also be a character string specifying a bandwidth selection method. Options available can be checked in \code{\link[stats:bandwidth]{bw.nrd}}. By default it is set to \verb{1/8 n^(-1/9)}, where \code{n} is the sample size.}

\item{lower}{Numeric value specifying the lower bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and sample size.}

\item{upper}{Numeric value specifying the upper bound for bandwidth selection. Default is computed based on the interquartile range (IQR) and sample size.}

\item{nh}{An integer specifying the number of points in the grid to evaluate the mean integrated squared error function. Default is 200.}

\item{tol}{Tolerance value used to check whether the minimum found lies at the boundaries of the interval; that is, the function will return a warning if the window minimizing the cross-validation function lies within \verb{[lower, lower+tol]} or \verb{[upper-tol, upper]}. Default is 10\% of the lower bound.}

\item{from}{Numeric value specifying the lower bound to be used in \code{\link[stats:density]{density}}. Default is computed based on the range of input data.}

\item{to}{Numeric value specifying the upper bound to be used in \code{\link[stats:density]{density}}. Default is computed based on the range of input data.}

\item{plot}{Logical value indicating whether to plot the mean integrated squared error function. Default is \code{TRUE}.}
}
\value{
The bootstrap bandwidth value.
}
\description{
This function computes the bandwidth selector for \insertCite{jones1991;textual}{WData} kernel density estimator using the bias-corrected bootstrap method developed by \insertCite{borrajo2017;textual}{WData}.
}
\details{
The bandwidth returned is the one minimizing \eqn{\mathrm{MISE}^{\ast}} over a compact interval \eqn{[h_1,h_2]} (determined by arguments \code{lower} and \code{upper}), i.e.,
\deqn{
\widehat{h}_{f, \mathrm{B}} = \arg \min_{h_{f} \in [h_1,h_2]} \int_{-\infty}^{+ \infty} \mathrm{MSE}^{\ast}\left(\widehat{f}^{\ast}_{\mathrm{J}, h_{f}}(y)\right)dy.
}
\eqn{\mathrm{MISE}^{\ast}} and \eqn{\mathrm{MSE}^{\ast}} correspond with the expression of the mean integrated squared error and the mean squared error of the bootstrap estimator \eqn{\widehat{f}^{\ast}_{\mathrm{J}, h_{f}}} provided by \insertCite{borrajo2017;textual}{WData}.
}
\examples{
\donttest{
bw.f.BGM.boot2(shrub.data$Width, nh = 50L)
}
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[=df.jones]{df.jones}}
}
