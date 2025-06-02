#' Generate a biased dataset using \insertCite{neumann1951;textual}{WData} acceptance-rejection method
#'
#' This function generates a biased dataset of size `n` using \insertCite{neumann1951;textual}{WData} acceptance-rejection method. The generated dataset is biased according to the provided bias function `w`, with respect to the unbiased density function `fx`.
#'
#' @param n Number of data points to generate.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `y`. By default, it is set to the length-biased function.
#' @param fx Unbiased density function. Distributions allowed are
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
#' Mixture of normal distributions ([`mixnorm`][KScorrect::dmixnorm()]) and
#' Mixture of gamma distributions ([`mgamma`][evmix::dmgamma()]) .
#' @param lim Lower and upper limits for the range where the bias is significant and, hence, where \eqn{c = \max_{y\in \mathbb{R}}w(y)/\mu_w > 0} must be searched.
#' @param plot Logical value indicating whether to generate a diagnostic plot of the biased dataset. Default is `TRUE`.
#' @param stop Logical value indicating whether to stop when bias function can't be evaluated in a generated value. Default is `TRUE`. If `FALSE` value is discarded and a new one is generated.
#' @param shape1,shape2 Additional arguments to be passed to the unbiased density function `fx` when set to the [`beta`][stats::dbeta()] distribution.
#' @param df,ncp Additional arguments to be passed to the unbiased density function `fx` when set to the [`chisq`][stats::dchisq()] and [`t`][stats::dt()] distributions.
#' @param shape,rate,scale,location Additional arguments to be passed to the unbiased density function `fx` when set to the [`cauchy`][stats::dcauchy()], [`logis`][stats::dlogis()], [`exp`][stats::dexp()], [`gamma`][stats::dgamma()] and [`weibull`][stats::dweibull()] distributions.
#' @param df1,df2 Additional arguments to be passed to the unbiased density function `fx` when set to the [`f`][stats::df()] distribution.
#' @param min,max Additional arguments to be passed to the unbiased density function `fx` when set to the [`unif`][stats::dunif()] distribution.
#' @param mean,sd,pro,meanlog,sdlog Additional arguments to be passed to the unbiased density function `fx` when set to the [`norm`][stats::dnorm()],  [`mixnorm`][KScorrect::dmixnorm()] and [`lnorm`][stats::dlnorm()]  distributions.
#' @param mgshape,mgscale,mgweight Additional arguments to be passed to the unbiased density function `fx` when set to the [`mgamma`][evmix::dmgamma()] distribution.
#' @return A numeric vector containing biased random samples from density `fx` and bias function `w`.
#' @details This function implements \insertCite{neumann1951;textual}{WData} acceptance-rejection method to generate a biased dataset given an  unbiased density function `fx` and a bias function `w`.
#' @references \insertAllCited{}
#' @examples
#' # Generate a length-biased dataset of size 100 from an exponential distribution
#' rbiased(n = 100, fx = "exp", rate = 2, plot = FALSE)
#'
#' #' # Generate a length-biased biased dataset from a gamma distribution
#' rbiased(n = 100, fx = "gamma", rate = 1.5^2, shape = 1.5)
#'
#' # Generate a biased dataset from a normal distribution
#' custom_bias <- function(y) {
#'   y^2
#' }
#' rbiased(n = 100, w = custom_bias, fx = "norm", mean = 3, sd = 10, plot = TRUE)
#'
#' # Generate a biased dataset from a mixture of normal distributions
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
    beta = function(x) {
      dbeta(x, shape1, shape2)
    },
    cauchy = function(x) {
      dcauchy(x, location, scale)
    },
    chisq = function(x) {
      dchisq(x, df, ncp)
    },
    exp = function(x) {
      dexp(x, rate)
    },
    f = function(x) {
      df(x, df1, df2)
    },
    gamma = function(x) {
      dgamma(x, shape, rate)
    },
    logis = function(x) {
      dlogis(x, location, scale)
    },
    lnorm = function(x) {
      dlnorm(x, meanlog, sdlog)
    },
    norm = function(x) {
      dnorm(x, mean, sd)
    },
    t = function(x) {
      dt(x, df, ncp)
    },
    unif = function(x) {
      dunif(x, min, max)
    },
    weibull = function(x) {
      dweibull(x, shape, scale)
    },
    mixnorm = function(x) {
      KScorrect::dmixnorm(x, mean, sd, pro)
    },
    mgamma = function(x) {
      evmix::dmgamma(x, mgshape, mgscale, mgweight)
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
    x <- rep(0, n)
    for (k in 1:G) {
      x[clabels == k] <- mean[k] + rnorm(ctable[k], sd = sd[k])
    }
    structure(as.vector(x), modelName = modelName)
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
    x <- rep(0, n)

    for (i in 1:G) {
      x[clabels == i] <- rgamma(1, shape = mgshape[i], scale = mgscale[i])
    }

    structure(as.vector(x), modelName = modelName)
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
          stop("bias function 'w' must be evaluable and positive in each point of 'fx' domain"),
          warning("bias function 'w' must be evaluable and positive in each point of 'fx' domain")
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
    stop("failure on the method -  it was not possible to generate biased dataset")
  }

  boundw <- sapply(seq(base::min(y, na.rm = T), base::max(y, na.rm = T), length.out = 1000L), w)
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
      w(y) * f(x = y)
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
      ylim = c(y_min, y_max),
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
