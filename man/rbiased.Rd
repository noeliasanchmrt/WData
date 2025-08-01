% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbiased.R
\name{rbiased}
\alias{rbiased}
\title{Generate a biased sample using \insertCite{neumann1951;textual}{WData} acceptance-rejection method}
\usage{
rbiased(
  n,
  w = function(y) {
     ifelse(y >= 0, y, NA)
 },
  fx,
  lim = 0.01,
  plot = TRUE,
  stop = TRUE,
  shape1,
  shape2,
  location,
  scale,
  df,
  ncp,
  rate,
  df1,
  df2,
  shape,
  meanlog,
  sdlog,
  min,
  max,
  mgshape,
  mgscale,
  mgweight,
  pro,
  mean,
  sd
)
}
\arguments{
\item{n}{Sample size.}

\item{w}{A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample \code{y}. By default, it is set to the length-biased function.}

\item{fx}{Unbiased density function. Values allowed are
Beta (\code{\link[stats:Beta]{beta}}),
Cauchy (\code{\link[stats:Cauchy]{cauchy}}),
Chi-Square (\code{\link[stats:Chisquare]{chisq}}),
Exponential (\code{\link[stats:Exponential]{exp}}),
F (\code{\link[stats:Fdist]{f}}),
Gamma  (\code{\link[stats:GammaDist]{gamma}}),
Logistic  (\code{\link[stats:Logistic]{logis}}),
Log Normal (\code{\link[stats:Lognormal]{lnorm}}),
Normal (\code{\link[stats:Normal]{norm}}),
Student t (\code{\link[stats:TDist]{t}}),
Continuous Uniform (\code{\link[stats:Uniform]{unif}}),
Weibull (\code{\link[stats:Weibull]{weibull}}),
Mixture of gaussian distributions (\code{\link[KScorrect:dmixnorm]{mixnorm}}) and
Mixture of gamma distributions (\code{\link[evmix:mgamma]{mgamma}}) .}

\item{lim}{Lower and upper limits for the range where the bias is significant and, hence, where \eqn{c = \max_{y\in \mathbb{R}}w(y)/\mu_w > 0} must be searched.}

\item{plot}{Logical value indicating whether to generate a plot of the biased sample. Default is \code{TRUE}.}

\item{stop}{Logical value indicating whether to stop when bias function can not be evaluated in a generated value. Default is \code{TRUE}. If \code{FALSE} value is discarded and a new one is generated.}

\item{shape1, shape2}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Beta]{beta}} distribution.}

\item{df, ncp}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Chisquare]{chisq}} and \code{\link[stats:TDist]{t}} distributions.}

\item{df1, df2}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Fdist]{f}} distribution.}

\item{shape, rate, scale, location}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Cauchy]{cauchy}}, \code{\link[stats:Logistic]{logis}}, \code{\link[stats:Exponential]{exp}}, \code{\link[stats:GammaDist]{gamma}} and \code{\link[stats:Weibull]{weibull}} distributions.}

\item{min, max}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Uniform]{unif}} distribution.}

\item{mgshape, mgscale, mgweight}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[evmix:mgamma]{mgamma}} distribution.}

\item{mean, sd, pro, meanlog, sdlog}{Additional arguments to be passed to the unbiased density function \code{fx} when set to the \code{\link[stats:Normal]{norm}},  \code{\link[KScorrect:dmixnorm]{mixnorm}} and \code{\link[stats:Lognormal]{lnorm}}  distributions.}
}
\value{
A numeric vector containing a biased sample from density \code{fx} and bias function \code{w}.
}
\description{
This function generates a biased sample of size \code{n} using \insertCite{neumann1951;textual}{WData} acceptance-rejection method. The generated sample is biased according to the provided bias function \code{w}, with respect to the unbiased density function \code{fx}.
}
\details{
This function implements \insertCite{neumann1951;textual}{WData} acceptance-rejection method to generate a biased sample given an  unbiased density function \code{fx} and a bias function \code{w}.
}
\examples{
# Generate a length-biased sample of size 100 from an exponential distribution
rbiased(n = 100, fx = "exp", rate = 2, plot = FALSE)

# Generate a length-biased sample from a gamma distribution
rbiased(n = 100, fx = "gamma", rate = 1.5^2, shape = 1.5)

# Generate a biased sample from a gaussian distribution
custom_bias <- function(y) {
  y^2
}
rbiased(n = 100, w = custom_bias, fx = "norm", mean = 3, sd = 10, plot = TRUE)

# Generate a biased sample from a mixture of gaussian distributions
custom_bias <- function(y) {
  sqrt(abs(y)) + 5
}
rbiased(
  n = 100, w = custom_bias, fx = "mixnorm", pro = rep(1 / 3, 3), mean = c(0.25, 0.5, 0.75),
  sd = rep(0.075, 3)
)
}
\references{
\insertAllCited{}
}
