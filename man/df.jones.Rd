% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df.jones.R
\name{df.jones}
\alias{df.jones}
\title{\insertCite{jones1991;textual}{WData} kernel density estimator}
\usage{
df.jones(
  y,
  w = function(y) {
     ifelse(y >= 0, y, NA)
 },
  y.seq,
  bw = "bw.f.BGM.rt",
  kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight",
    "cosine", "optcosine"),
  from,
  to,
  nb = 512L,
  plot = TRUE,
  ...
)
}
\arguments{
\item{y}{A numeric vector containing the biased sample.}

\item{w}{A function representing the bias function applied to the data points.
It must be evaluable and positive in each point of the sample \code{y}.
By default, it is set to the length-biased function.}

\item{y.seq}{A numeric vector specifying the points where the density is estimated. Alternatively, \code{from}, \code{to} and \code{nb} can be used to define the evaluation points.}

\item{bw}{The smoothing bandwidth to be used in the density estimation. \code{bw} can also be a character string giving a rule to choose the bandwidth. In this case, options available are \code{\link[=bw.f.BGM.rt]{bw.f.BGM.rt}}, \code{\link[=bw.f.BGM.cv]{bw.f.BGM.cv}},  \code{\link[=bw.f.BGM.boot1]{bw.f.BGM.boot1}} and \code{\link[=bw.f.BGM.boot2]{bw.f.BGM.boot2}}.
Default is \code{\link[=bw.f.BGM.rt]{bw.f.BGM.rt}}.}

\item{kernel}{A character string specifying the kernel function. Available options: \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, \code{"biweight"}, \code{"cosine"} and \code{"optcosine"}.}

\item{from}{Numeric value specifying the lower bound of the grid where the estimator is computed when \code{y.seq} is not provided. Default is computed based on the range of input data.}

\item{to}{Numeric value specifying the upper bound of the grid where the estimator is computed when \code{y.seq} is not provided. Default is computed based on the range of input data.}

\item{nb}{An integer specifying the number of points at which the estimator is computed when \code{y.seq} is not provided. Default is 512.}

\item{plot}{A logical value indicating whether to plot the density estimation. Default is \code{TRUE}.}

\item{...}{Additional arguments to be passed to bandwidth selection functions.}
}
\value{
A list with the following components:
\item{\code{y.seq}}{The points where the density is estimated.}
\item{\code{f.hat}}{The estimated density values.}
\item{\code{bw}}{The bandwidth value.}
\item{\code{n}}{The sample size after removal of \code{NaN}, \code{Na} and \code{Inf}.}
\item{\code{call}}{The call which produced the result.}
\item{\code{has.na}}{Logical; indicates whether the original vector \code{y} contains any \code{NaN}, \code{Na} or \code{Inf}.}
}
\description{
This function computes \insertCite{jones1991;textual}{WData} kernel density estimator given a sample and the corresponding biased function.
}
\details{
\insertCite{jones1991;textual}{WData} kernel density estimator is expressed as
\deqn{\widehat{f}_{\mathrm{J}, h_{f}}(y)=\frac{\widehat{\mu}_w}{n}\sum_{i=1}^{n}  \frac{1}{w(Y_i)} K_{h_{f}}(y-Y_i),\quad
\text{where}\quad \widehat{\mu}_w=n \left(\sum_{i=1}^{n} \frac{1}{w(Y_i)}\right)^{-1},}
\eqn{h_{f}} is the bandwidth, \eqn{K} is the kernel density function and \eqn{K_{h_{f}}(u)=1/h_{f} K\left(u / h_{f}\right)}.
}
\examples{
# Rule of thumb
df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.rt")
# Cross Validation
df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.cv")
# Bootstrap
df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot1", bw0 = "RT")
df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot1", bw0 = "PI")
\donttest{
df.jones(y = shrub.data$Width, kernel = "epanechnikov", bw = "bw.f.BGM.boot2", nh = 50L)
}
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[=bw.f.BGM.rt]{bw.f.BGM.rt}}, \code{\link[=bw.f.BGM.cv]{bw.f.BGM.cv}},  \code{\link[=bw.f.BGM.boot1]{bw.f.BGM.boot1}} , \code{\link[=bw.f.BGM.boot2]{bw.f.BGM.boot2}}
}
