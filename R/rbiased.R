#' Generate a biased sample using \insertCite{neumann1951;textual}{WData} acceptance-rejection method
#'
#' This function generates a biased sample of size `n` using \insertCite{neumann1951;textual}{WData} acceptance-rejection method. The generated sample is biased according to the provided bias function `w`, with respect to the unbiased density function `fx`.
#'
#' @param n Sample size.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of the sample `y`. By default, it is set to the length-biased function.
#' @param fx Unbiased density function. Values allowed are
#' Beta ([`beta`][stats::dbeta()]),
#' Cauchy ([`cauchy`][stats::dcauchy()]),
#' Chi-Square ([`chisq`][stats::dchisq()]),
#' Exponential ([`exp`][stats::dexp()]),
#' F ([`f`][stats::df()]),
#' Gamma  ([`gamma`][stats::dgamma()]),
#' Logistic  ([`logis`][stats::dlogis()]),
#' Log Normal ([`lnorm`][stats::dlnorm()]),
#' Normal ([`norm`][stats::dnorm()]),
#' Student t ([`t`][stats::dt()]),
#' Continuous Uniform ([`unif`][stats::dunif()]),
#' Weibull ([`weibull`][stats::dweibull()]),
#' Mixture of gaussian distributions ([`mixnorm`][KScorrect::dmixnorm()]) and
#' Mixture of gamma distributions ([`mgamma`][evmix::dmgamma()]) .
#' @param lim A numeric value between 0 and 1 specifying the quantile range within which the bias function is considered significant and, hence, the interval over which the constant \eqn{c = \max_{y \in \mathbb{R}} w(y)/\mu_w > 0} is searched. The lower and upper bounds for this interval are determined by `q(lim)` and `q(1 - lim)``, respectively, where `q` is the quantile function of the chosen distribution `fx`.
#' @param lim Lower and upper limits for the range where the bias is significant and, hence, where \eqn{c = \max_{y\in \mathbb{R}}w(y)/\mu_w > 0} must be searched.
#' @param plot Logical value indicating whether to generate a plot of the biased sample. Default is `TRUE`.
#' @param stop Logical value indicating whether to stop when bias function can not be evaluated in a generated value. Default is `TRUE`. If `FALSE` value is discarded and a new one is generated.
#' @param shape1,shape2 Additional arguments to be passed to the unbiased density function `fx` when set to the [`beta`][stats::dbeta()] distribution.
#' @param df,ncp Additional arguments to be passed to the unbiased density function `fx` when set to the [`chisq`][stats::dchisq()] and [`t`][stats::dt()] distributions.
#' @param shape,rate,scale,location Additional arguments to be passed to the unbiased density function `fx` when set to the [`cauchy`][stats::dcauchy()], [`logis`][stats::dlogis()], [`exp`][stats::dexp()], [`gamma`][stats::dgamma()] and [`weibull`][stats::dweibull()] distributions.
#' @param df1,df2 Additional arguments to be passed to the unbiased density function `fx` when set to the [`f`][stats::df()] distribution.
#' @param min,max Additional arguments to be passed to the unbiased density function `fx` when set to the [`unif`][stats::dunif()] distribution.
#' @param mean,sd,pro,meanlog,sdlog Additional arguments to be passed to the unbiased density function `fx` when set to the [`norm`][stats::dnorm()],  [`mixnorm`][KScorrect::dmixnorm()] and [`lnorm`][stats::dlnorm()]  distributions.
#' @param mgshape,mgscale,mgweight Additional arguments to be passed to the unbiased density function `fx` when set to the [`mgamma`][evmix::dmgamma()] distribution.
#' @return A numeric vector containing a biased sample from density `fx` and bias function `w`.
#' @details This function implements \insertCite{neumann1951;textual}{WData} acceptance-rejection method to generate a biased sample given an  unbiased density function `fx` and a bias function `w`.
#' @references \insertAllCited{}
#' @export
#' @examples
#' # Generate a length-biased sample of size 100 from an exponential distribution
#' rbiased(n = 100, fx = "exp", rate = 2, plot = FALSE)
#'
#' # Generate a length-biased sample from a gamma distribution
#' rbiased(n = 100, fx = "gamma", rate = 1.5^2, shape = 1.5)
#'
#' # Generate a biased sample from a gaussian distribution
#' custom_bias <- function(y) {
#'   y^2
#' }
#' rbiased(n = 100, w = custom_bias, fx = "norm", mean = 3, sd = 10, plot = TRUE)
#'
#' # Generate a biased sample from a mixture of gaussian distributions
#' custom_bias <- function(y) {
#'   sqrt(abs(y)) + 5
#' }
#' rbiased(
#'   n = 100, w = custom_bias, fx = "mixnorm", pro = rep(1 / 3, 3), mean = c(0.25, 0.5, 0.75),
#'   sd = rep(0.075, 3)
#' )
rbiased <- function(n,
                    w = function(y) {
                      ifelse(y >= 0, y, NA)
                    },
                    fx,
                    lim = 0.01,
                    plot = TRUE,
                    stop = TRUE,
                    shape1, shape2,
                    location, scale,
                    df, ncp, rate, df1, df2,
                    shape, meanlog, sdlog, min, max,
                    mgshape, mgscale, mgweight, pro, mean, sd) {
  # check_dots_used()
  if (!isTRUE(all.equal(n, as.integer(n))) | n < 1) {
    stop("argument 'n' must be and integer greater or equal than 1")
  }
  if (!is.function(w)) {
    stop("argument 'w' must be a function")
  }
  if (!is.logical(plot)) {
    stop("argument 'plot' must be logical")
  }
  if (!is.numeric(lim) | length(lim) != 1 | lim <= 0 | lim >= 1) {
    stop("'lim' must be an numeric value greater than 0 and small than 1")
  }

  fx <- match.arg(fx, c(
    "beta", "cauchy", "chisq", "exp", "f", "gamma", "logis", "lnorm", "norm", "t",
    "unif", "weibull", "mixnorm", "mgamma"
  ))

  f <- switch(fx,
    beta = function(y) {
      dbeta(y, shape1, shape2)
    },
    cauchy = function(y) {
      dcauchy(y, location, scale)
    },
    chisq = function(y) {
      dchisq(y, df, ncp)
    },
    exp = function(y) {
      dexp(y, rate)
    },
    f = function(y) {
      df(y, df1, df2)
    },
    gamma = function(y) {
      dgamma(y, shape, rate)
    },
    logis = function(y) {
      dlogis(y, location, scale)
    },
    lnorm = function(y) {
      dlnorm(y, meanlog, sdlog)
    },
    norm = function(y) {
      dnorm(y, mean, sd)
    },
    t = function(y) {
      dt(y, df, ncp)
    },
    unif = function(y) {
      dunif(y, min, max)
    },
    weibull = function(y) {
      dweibull(y, shape, scale)
    },
    mixnorm = function(y) {
      KScorrect::dmixnorm(y, mean, sd, pro)
    },
    mgamma = function(y) {
      evmix::dmgamma(y, mgshape, mgscale, mgweight)
    }
  )

  q <- switch(fx,
    beta = function(p) {
      qbeta(p, shape1, shape2)
    },
    cauchy = function(p) {
      qcauchy(p, location, scale)
    },
    chisq = function(p) {
      qchisq(p, df, ncp)
    },
    exp = function(p) {
      qexp(p, rate)
    },
    f = function(p) {
      qf(p, df1, df2)
    },
    gamma = function(p) {
      qgamma(p, shape, rate)
    },
    logis = function(p) {
      qlogis(p, location, scale)
    },
    lnorm = function(p) {
      qlnorm(p, meanlog, sdlog)
    },
    norm = function(p) {
      qnorm(p, mean, sd)
    },
    t = function(p) {
      qt(p, df, ncp)
    },
    unif = function(p) {
      qunif(p, min, max)
    },
    weibull = function(p) {
      qweibull(p, shape, scale)
    },
    mixnorm = function(p) {
      KScorrect::qmixnorm(p, mean, sd, pro)
    },
    mgamma = function(p) {
      evmix::qmgamma(p, mgshape, mgscale, mgweight)
    }
  )

  # Auxiliar Function to Generate a Mixture of Normals

  r_mixto_norm <- function(n, mean, sd, pro) {
    if (any(missing(mean), missing(sd))) stop("'mean' and 'sd' not provided, without default")
    mean <- as.vector(mean, mode = "numeric")
    G <- length(mean)
    sd <- as.vector(sd, mode = "numeric")

    if (missing(pro)) {
      pro <- rep(1 / G, G)
      warning("mixing proportion 'pro' not provided. Assigned equal proportions by default")
    }
    if (any(pro < 0L, sd < 0L)) stop("'pro' and 'sd' must not be negative")
    lpro <- length(pro)
    modelName <- "V"
    lsd <- length(sd)
    if (lsd == 1L & G > 1L) {
      modelName <- "E"
      sd[seq(G)] <- sd[1]
      lsd <- length(sd)
      warning("'equal variance model' implemented. If want 'variable-variance model', specify remaining 'sd's")
    }
    if (G < lsd | G < lpro | (lsd > 1L & G != lsd) | (!missing(pro) & G != lpro)) {
      stop("the lengths of supplied parameters do not make sense")
    }
    pro <- as.vector(pro, mode = "numeric")
    pro <- pro / sum(pro)
    clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
    ctable <- tabulate(clabels, nbins = G)
    y <- rep(0, n)
    for (k in 1:G) {
      y[clabels == k] <- mean[k] + rnorm(ctable[k], sd = sd[k])
    }
    structure(as.vector(y), modelName = modelName)
  }


  # Auxiliar Function to Generate a Mixture of gammas
  r_mixto_gamma <- function(n, mgshape, mgscale, mgweight) {
    if (any(missing(mgshape), missing(mgscale))) {
      stop("'mgshape' and 'mgscale' not provided, without default")
    }

    if (length(mgshape) != length(mgweight)) {
      stop("lengths of 'shape', 'scale', and 'weight' must be equal")
    }
    G <- length(mgshape)
    if (missing(mgweight)) {
      mgweight <- rep(1 / G, G)
      warning("mixing proportion 'mgweight' not provided. Assigned equal proportions by default")
    }

    mgweight <- as.vector(mgweight, mode = "numeric")
    mgweight <- mgweight / sum(mgweight)
    if (any(mgshape < 0L, mgscale < 0L, mgweight < 0L)) {
      stop("'mgshape', 'mgscale' and 'mgweight' must not be negative")
    }

    modelName <- "V"
    clabels <- sample(1:G, size = n, replace = TRUE, prob = mgweight)
    ctable <- tabulate(clabels, nbins = G)
    y <- rep(0, n)

    for (i in 1:G) {
      y[clabels == i] <- rgamma(1, shape = mgshape[i], scale = mgscale[i])
    }

    structure(as.vector(y), modelName = modelName)
  }

  r <- switch(fx,
    beta = function(n) {
      rbeta(n, shape1, shape2)
    },
    cauchy = function(n) {
      rcauchy(n, location, scale)
    },
    chisq = function(n) {
      rchisq(n, df, ncp)
    },
    exp = function(n) {
      rexp(n, rate = rate)
    },
    f = function(n) {
      rf(n, df1, df2)
    },
    gamma = function(n) {
      rgamma(n, shape, rate)
    },
    logis = function(n) {
      rlogis(n, location, scale)
    },
    lnorm = function(n) {
      rlnorm(n, meanlog, sdlog)
    },
    norm = function(n) {
      rnorm(n, mean, sd)
    },
    t = function(n) {
      rt(n, df, ncp)
    },
    unif = function(n) {
      runif(n, min, max)
    },
    weibull = function(n) {
      rweibull(n, shape, scale)
    },
    mixnorm = function(n) {
      r_mixto_norm(n, mean, sd, pro)
    },
    mgamma = function(n) {
      r_mixto_gamma(n, mgshape, mgscale, mgweight)
    }
  )

  c <- optimize(w, interval = c(q(p = lim), q(p = (1 - lim))), maximum = TRUE)$maximum

  if (!is.numeric(c) | abs(c) <= 0) {
    stop("failure on the method -  it was wot possible to find suitable 'c'")
  }

  y <- numeric(n)
  ngen <- 0
  ar1 <- function(c, w) {
    while (TRUE) {
      xi <- r(1)
      wxi <- w(xi)
      if (any(length(wxi) == 0 | !is.numeric(wxi) | wxi <= 0)) {
        ifelse(stop,
          stop("bias function 'w' must be evaluable and positive in each point of the sample 'fx' domain"),
          warning("bias function 'w' must be evaluable and positive in each point of the sample 'fx' domain")
        )
      }

      ngen <<- ngen + 1
      if (ifelse(length(xi) == 1 & all(is.numeric(wxi)) & all(wxi > 0), abs(c) * runif(1) <= wxi, FALSE)) {
        return(xi)
      }
    }
  }

  for (i in 1:n) {
    y[i] <- ar1(c = c, w = w)
  }


  if (any(length(y) == 0 | !is.numeric(y) | anyNA(y))) {
    stop("failure on the method -  it was not possible to generate biased sample")
  }

  boundw <- sapply(seq(base::min(y, na.rm = TRUE), base::max(y, na.rm = TRUE), length.out = 1000L), w)
  if (any(length(boundw) == 0) | anyNA(boundw)) {
    stop("function 'w' must be bounded in [min(y), max(y)]")
  }

  integrabilityw <- try(integrate(w, lower = base::min(y), upper = base::max(y))$value, silent = TRUE)
  if (inherits(integrabilityw, "try-error") &&
    grepl("the integral is probably divergent$", attr(integrabilityw, "condition")$message)) {
    stop("function 'w' must be integrable in [min(y), max(y)]")
  }

  if (plot == TRUE) {
    wf <- function(y) {
      w(y) * f(y = y)
    }
    uw <- integrate(wf, lower = q(lim), upper = q(1 - lim))
    g <- function(y) {
      wf(y) / uw$value
    }

    y_values <- seq(base::min(y), base::max(y), length.out = 1000)
    density_vals <- density(y, bw = bw.SJ(y))$y
    g_vals <- g(y_values)
    f_vals <- f(y_values)

    y_min <- base::min(c(density_vals, g_vals, f_vals), na.rm = TRUE)
    y_max <- base::max(c(density_vals, g_vals, f_vals), na.rm = TRUE)

    plot(density(y, bw = bw.SJ(y)),
      col = "black", lty = 2,
      ylim = c(y_min, 1.2 * y_max),
      sub = paste0("Bias:", gsub("\\s+", " ", deparse1(w))),
      xlab = paste0("Unbiased Density f(y):", gsub("\\s+", " ", deparse1(fx)))
    )
    rug(y, col = "black")
    curve(g, add = TRUE, col = "black")
    curve(f, add = TRUE, col = "blue")

    legend("topright",
      legend = c(
        "Simulated Biased Density",
        "Expected Simulated Biased Density",
        "Unbiased Density"
      ),
      col = c("black", "black", "blue"),
      lty = c(2, 1, 1, 1), cex = 0.8
    )
  }
  return(y)
}
