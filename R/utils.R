.check_biased_sample <- function(y, w) {
  data.name <- deparse(substitute(y))

  if (!is.numeric(y)) {
    stop("invalid 'y'")
  }

  has.na <- any(is.na(y))
  if ((n <- length(y <- y[!is.na(y)])) < 2L) {
    stop("need at least 2 data points")
  }

  if (is.na(n <- as.integer(n))) {
    stop("invalid length(y)")
  }

  if (!is.function(w)) {
    stop("argument 'w' must be a function")
  }

  boundw <- sapply(seq(min(y), max(y), length.out = 1000L), w)
  if (any(length(boundw) == 0) | anyNA(boundw)) {
    stop("function 'w' must be bounded in [min(y), max(y)]")
  }

  integrabilityw <- try(integrate(w, lower = min(y), upper = max(y))$value, silent = TRUE)
  if (inherits(integrabilityw, "try-error") &&
    grepl("the integral is probably divergent$", attr(integrabilityw, "condition")$message)) {
    stop("function 'w' must be integrable in [min(y), max(y)]")
  }

  yw <- sapply(y, w)
  if (any(length(yw) == 0) | anyNA(yw) | any(yw <= 0)) {
    stop("function 'w' must be evaluable and positive in each observation")
  }

  weights <- yw^(-1)

  uw <- sum(sapply(yw, function(t) {
    1 / t
  }))^(-1) * n
  if (!is.numeric(uw)) stop("non numeric value for 'uw'")

  uwb <- sum(sapply(yw, function(t) {
    1 / t
  }))^(-1) *
    sum(sapply(yw, function(t) {
      1 / t^2
    }))
  if (!is.numeric(uwb)) stop("non numeric value for 'uwb'")

  list(data.name = data.name, n = n, y = y, w = w, yw = yw, weights = weights, uw = uw, uwb = uwb, has.na = has.na)
}

.get_kernel_values <- function(kernel) {
  kernel_function_density <- # Kernel Density functions
    switch(kernel,
      gaussian = function(u) dnorm(u),
      epanechnikov = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, 0.75 * (1 - (u^2)))),
      rectangular = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, 0.5)),
      triangular = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, 1 - abs(u))),
      biweight = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, 15 / 16 * (1 - u^2)^2)),
      cosine = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, (1 + cos(pi * u)) / 2)),
      optcosine = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 0, pi / 4 * cos(pi * u / 2))),
      stop("unknown kernel")
    )

  kernel_function_distribution <- # Kernel Distribution functions
    switch(kernel,
      gaussian = function(u) pnorm(u),
      epanechnikov = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, 0.75 * u * (1 - (u^2) / 3) + 0.5)),
      rectangular = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, (u + 1) / 2)),
      triangular = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, -0.5 * u * abs(u) + u + 0.5)),
      biweight = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, (15 / 16) * u - (5 / 8) * u^3 + (3 / 16) * u^5 + 0.5)),
      cosine = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, (sin(pi * u) + pi * u + pi) / (2 * pi))),
      optcosine = function(u) ifelse(u <= -1, 0, ifelse(u >= 1, 1, (1 + sin(pi * u / 2)) / 2)),
      stop("unknown kernel")
    )

  kernel_function_density_deriv1 <- switch(kernel, # Kernel First Derivative of density functions
    gaussian = function(u) dnorm(u) * (-u),
    epanechnikov = function(u) ifelse(abs(u) < 1, -1.5 * u, 0),
    rectangular = function(u) ifelse(abs(u) < 1, 0, 0),
    triangular = function(u) ifelse(abs(u) < 1 & abs(u) > .Machine$double.eps^0.5, -u / abs(u), 0),
    biweight = function(u) ifelse(abs(u) < 1, -(15 / 4) * u * (1 - u^2), 0),
    cosine = function(u) -0.5 * sin(pi * u) * pi * (abs(u) <= 1),
    optcosine = function(u) -pi^2 / 8 * sin(pi / 2 * u) * (abs(u) <= 1),
    stop("unknown kernel")
  )

  kernel_function_density_deriv2 <- switch(kernel,
    gaussian = function(u) {
      1 / (sqrt(2 * pi)) * (u^2 - 1) * exp(-u^2 / 2)
    },
    epanechnikov = function(u) {
      -3 / 2 * (abs(u) <= 1)
    },
    rectangular = function(u) {
      0 * (abs(u) <= 1)
    },
    triangular = function(u) {
      0 * (abs(u) <= 1)
    },
    biweight = function(u) {
      0.25 * (45 * u^2 - 15) * (abs(u) <= 1)
    },
    cosine = function(u) {
      -0.5 * pi^2 * cos(pi * u) * (abs(u) <= 1)
    },
    optcosine = function(u) {
      -pi^3 / 16 * cos(pi * u / 2) * (abs(u) <= 1)
    },
    stop("unknown kernel")
  )

  if (kernel != "gaussian" & kernel != "rectangular") {
    eps <- switch(kernel,
      gaussian = NA,
      epanechnikov = 0.001,
      rectangular = NA,
      triangular = 0.0025,
      biweight = 0.001,
      cosine = 0.001,
      optcosine = 0.005,
      stop("unknown kernel")
    )

    conv <- bayesmeta::convolve(dens1 = kernel_function_density, dens2 = kernel_function_density, epsilon = eps)$density
    u_vals <- seq(-2, 2, length.out = 1000)
    interpolated_conv <- approxfun(u_vals, conv(u_vals), rule = 2)

    aux <- function(u) {
      result <- interpolated_conv(as.vector(u))
      matrix(result, nrow = nrow(u), ncol = ncol(u))
    }
  }

  kernel_function_conv <-
    switch(kernel,
      gaussian = function(u) dnorm(u, sd = sqrt(2)),
      epanechnikov = aux,
      rectangular = function(u) ifelse(abs(u) <= 2, 0.25 * pmin(1, u + 1) - 0.25 * pmax(-1, u - 1), 0),
      triangular = aux,
      biweight = aux,
      cosine = aux,
      optcosine = aux,
      stop("unknown kernel")
    )

  # kernel_r <- switch(kernel,
  #   gaussian = 1 / (2 * sqrt(pi)),
  #   epanechnikov = 3 / 5,
  #   rectangular = 1 / 2,
  #   triangular = 2 / 3,
  #   biweight = 5 / 7,
  #   cosine = 3 / 4,
  #   optcosine = pi^2 / 16,
  #   stop("unknown kernel")
  # )

  kernel_r <- integrate(
    Vectorize(function(u) kernel_function_density(u)^2),
    -4, 4
  )$value

  # kernel_r_deriv1 <- switch(kernel,
  #   gaussian = 1 / (4 * sqrt(pi)),
  #   epanechnikov = 3 / 2,
  #   rectangular = 0,
  #   triangular = 2,
  #   biweight = 15 / 7,
  #   cosine = pi^2 / 4,
  #   optcosine = pi^4 / 64,
  #   stop("unknown kernel")
  # )

  kernel_r_deriv1 <- integrate(
    Vectorize(function(u) kernel_function_density_deriv1(u)^2),
    -4, 4
  )$value

  # kernel_r_deriv2 <- switch(kernel,
  #   gaussian = 3 / (8 * sqrt(pi)),
  #   epanechnikov = 9 / 2,
  #   rectangular = 0,
  #   triangular = 0,
  #   biweight = 22.5,
  #   cosine = pi^4 / 4,
  #   optcosine = pi^6 / 256,
  #   stop("unknown kernel")
  # )

  kernel_r_deriv2 <- integrate(
    Vectorize(function(u) kernel_function_density_deriv2(u)^2),
    -4, 4
  )$value

  # kernel_eta <- switch(kernel,
  #   gaussian = 1,
  #   epanechnikov = 1 / 5,
  #   rectangular = 1 / 3,
  #   triangular = 1 / 6,
  #   biweight = 1 / 7,
  #   cosine = (pi^2 - 6) / (3 * pi^2),
  #   optcosine = (pi^2 - 8) / (pi^2),
  #   stop("unknown kernel")
  # )

  kernel_eta <- integrate(
    Vectorize(function(u) u^2 * kernel_function_density(u)),
    -4, 4
  )$value

  # kernel_kappa <- switch(kernel,
  #   gaussian = 1 / sqrt(pi),
  #   epanechnikov = 9 / 35,
  #   rectangular = 1 / 3,
  #   triangular = 7 / 30,
  #   biweight = 50 / 231,
  #   cosine = (4 * pi^2 - 15) / (12 * pi^2),
  #   optcosine = 1 / 4,
  #   stop("unknown kernel")
  # )

  kernel_kappa <- integrate(
    Vectorize(function(u) 2 * u * kernel_function_density(u) * kernel_function_distribution(u)),
    -4, 4
  )$value

  list(
    kernel_function_density = kernel_function_density,
    kernel_function_distribution = kernel_function_distribution,
    kernel_function_density_deriv1 = kernel_function_density_deriv1,
    kernel_function_density_deriv2 = kernel_function_density_deriv2,
    kernel_function_conv = kernel_function_conv,
    kernel_r = kernel_r,
    kernel_r_deriv1 = kernel_r_deriv1,
    kernel_r_deriv2 = kernel_r_deriv2,
    kernel_eta = kernel_eta,
    kernel_kappa = kernel_kappa
  )
}


.simpsons_rule <- function(y.seq, fx) {
  # https://scicomp.stackexchange.com/questions/25649/composite-simpsons-rule-with-odd-intervals
  valid <- which(!is.na(y.seq) & !is.na(fx) & is.finite(y.seq) & is.finite(fx))
  y.seq <- y.seq[valid]
  fx <- fx[valid]

  ord <- order(y.seq)
  y.seq <- y.seq[ord]
  fx <- fx[ord]
  n <- length(y.seq)
  if (n < 5) stop("At least 5 points are required for Simpson's rule")

  h <- (y.seq[2] - y.seq[1])
  if (any(abs(diff(y.seq) - h) > .Machine$double.eps^0.5)) stop("y.seq must be equally spaced")

  integral <- fx[1] + fx[n] + 4 * sum(fx[seq.int(2, n - 1, by = 2)]) + 2 * sum(fx[seq.int(3, n - 2, by = 2)])
  integral * h / 3
}


.validate_number <- function(x, name) {
  if (!is.numeric(x) || any(is.na(x)) || any(!is.finite(x)) || length(x) != 1) {
    stop(paste("Invalid", name, ":", name, "must be numeric of length one"))
  }
}

.get_bandwidth_grid <- function(nh, lower, upper, tol, plot) {
  if (!all.equal(nh, as.integer(nh)) || nh <= 0) stop("'nh' must be a positive integer")
  nh <- as.integer(nh)
  .validate_number(nh, "nh")
  .validate_number(lower, "lower")
  .validate_number(upper, "upper")
  .validate_number(tol, "tol")

  if (lower < 0 || lower < 0) stop("'from' and 'to' must be positive numbers")
  if (lower > upper) stop("Invalid interval: 'lower' must be smaller than 'upper'")
  if (lower + tol >= upper - tol) stop("Invalid interval: 'lower + tol' must be smaller than 'upper - tol'")
  if (!is.logical(plot) || length(plot) != 1) stop("Invalid 'plot': must be a single logical value (TRUE or FALSE)")

  message(sprintf("Interval where bandwidth is searched: [%f, %f]", lower, upper))
  seq(lower, upper, length.out = nh)
}

.get_xaxn_grid <- function(y, y.seq, from, to, nb, plot) {
  if (missing(y.seq)) {
    if (missing(from)) {
      from <- min(y) - (sort(y)[5] - min(y))
    }
    if (missing(to)) {
      to <- max(y) + (max(y) - sort(y, decreasing = TRUE)[5])
    }

    .validate_number(nb, "nb")
    .validate_number(from, "from")
    .validate_number(to, "to")
    if (from >= to) stop("'from' must be smaller than 'to'")
    if (!all.equal(nb, as.integer(nb)) || nb <= 0) stop("'nb' must be a positive integer")
    y.seq <- seq.int(from, to, length.out = nb)
  } else {
    if (!is.numeric(y.seq) || !is.vector(y.seq)) stop("'y.seq' must be a numeric vector")
    from <- min(y.seq)
    to <- max(y.seq)
  }

  message(sprintf("Interval for Estimation: [%f, %f]", min(y.seq), max(y.seq)))

  if (!is.logical(plot)) {
    stop("argument 'plot' must be logical")
  }
  list(y.seq = y.seq, from = from, to = to)
}

# Adapted from evmix 2.12 (GPL-3).
.gamma_mix_parameters <- function(shape, scale, weight = NULL) {
  if (missing(shape)) stop("'mgshape' not provided, without default")
  if (missing(scale)) stop("'mgscale' not provided, without default")

  shape <- as.numeric(shape)
  scale <- as.numeric(scale)
  if (is.null(weight)) {
    weight <- rep(1, length(shape))
  }
  weight <- as.numeric(weight)

  if (length(shape) == 0L || length(shape) != length(scale) || length(shape) != length(weight)) {
    stop("lengths of 'shape', 'scale', and 'weight' must be equal")
  }
  if (any(!is.finite(c(shape, scale, weight))) || any(shape <= 0) ||
    any(scale <= 0) || any(weight <= 0)) {
    stop("'shape', 'scale', and 'weight' must be positive and finite")
  }

  list(shape = shape, scale = scale, weight = weight / sum(weight))
}

# Adapted from evmix 2.12 (GPL-3).
.dgamma_mix <- function(x, shape, scale, weight = NULL, log = FALSE) {
  mix <- .gamma_mix_parameters(shape, scale, weight)
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("'log' must be a single logical value")
  }
  if (any(is.infinite(x))) warning("infinite quantiles set to NA")
  x[is.infinite(x)] <- NA_real_

  density <- x
  finite <- which(is.finite(x))
  if (length(finite) > 0L) {
    density[finite] <- 0
    for (i in seq_along(mix$shape)) {
      density[finite] <- density[finite] +
        stats::dgamma(x[finite], shape = mix$shape[i], scale = mix$scale[i]) * mix$weight[i]
    }
  }
  if (log) density <- log(density)
  density
}

# Adapted from evmix 2.12 (GPL-3).
.pgamma_mix <- function(q, shape, scale, weight = NULL, lower.tail = TRUE) {
  mix <- .gamma_mix_parameters(shape, scale, weight)
  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value")
  }
  if (any(is.infinite(q))) warning("infinite quantiles set to NA")
  q[is.infinite(q)] <- NA_real_

  probability <- q
  finite <- which(is.finite(q))
  if (length(finite) > 0L) {
    probability[finite] <- 0
    for (i in seq_along(mix$shape)) {
      probability[finite] <- probability[finite] +
        stats::pgamma(q[finite], shape = mix$shape[i], scale = mix$scale[i]) * mix$weight[i]
    }
  }
  if (!lower.tail) probability <- 1 - probability
  probability
}

# Adapted from evmix 2.12 (GPL-3).
.qgamma_mix <- function(p, shape, scale, weight = NULL, lower.tail = TRUE) {
  mix <- .gamma_mix_parameters(shape, scale, weight)
  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.numeric(p) || length(p) == 0L || any(is.infinite(p)) ||
    any(p < 0, na.rm = TRUE) || any(p > 1, na.rm = TRUE)) {
    stop("'p' must contain probabilities between 0 and 1")
  }
  if (!lower.tail) p <- 1 - p

  quantile <- p
  finite <- which(is.finite(p))
  if (length(finite) > 0L) {
    maxp <- max(ifelse(p[finite] >= 1, NA_real_, p[finite]), na.rm = TRUE)
    maxq <- max(stats::qgamma(maxp, shape = mix$shape, scale = mix$scale))
    qgrid <- seq(0, maxq, length.out = 1000L)
    pgrid <- .pgamma_mix(qgrid, mix$shape, mix$scale, mix$weight)
    qfun <- stats::splinefun(x = pgrid, y = qgrid)
    quantile[finite] <- qfun(p[finite])
  }
  quantile
}

# Adapted from evmix 2.12 (GPL-3)
.rgamma_mix <- function(n, shape, scale, weight = NULL) {
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || !is.finite(n) ||
    n < 1 || abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    stop("'n' must be a positive integer")
  }
  mix <- .gamma_mix_parameters(shape, scale, weight)
  component <- sample.int(length(mix$shape), size = n, replace = TRUE, prob = mix$weight)
  sample <- numeric(n)
  for (i in seq_along(mix$shape)) {
    selected <- which(component == i)
    sample[selected] <- stats::rgamma(
      length(selected),
      shape = mix$shape[i], scale = mix$scale[i]
    )
  }
  sample
}
