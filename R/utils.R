.check_biased_dataset <- function(y, w) {
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

  kernel_function_density_deriv <- switch(kernel, # Kernel First Derivative of density functions
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

  # RK <- switch(kernel,
  #   gaussian = 1 / (2 * sqrt(pi)),
  #   epanechnikov = 3 / 5,
  #   rectangular = 1 / 2,
  #   triangular = 2 / 3,
  #   biweight = 5 / 7,
  #   cosine = 3 / 4,
  #   optcosine = pi^2 / 16,
  #   stop("unknown kernel")
  # )

  RK <- integrate(
    Vectorize(function(u) kernel_function_density(u)^2),
    -4, 4
  )$value

  # RKprime <- switch(kernel,
  #   gaussian = 1 / (4 * sqrt(pi)),
  #   epanechnikov = 3 / 2,
  #   rectangular = 0,
  #   triangular = 2,
  #   biweight = 15 / 7,
  #   cosine = pi^2 / 4,
  #   optcosine = pi^4 / 64,
  #   stop("unknown kernel")
  # )

  RKprime <- integrate(
    Vectorize(function(u) kernel_function_density_deriv(u)^2),
    -4, 4
  )$value

  # RKprime2 <- switch(kernel,
  #   gaussian = 3 / (8 * sqrt(pi)),
  #   epanechnikov = 9 / 2,
  #   rectangular = 0,
  #   triangular = 0,
  #   biweight = 22.5,
  #   cosine = pi^4 / 4,
  #   optcosine = pi^6 / 256,
  #   stop("unknown kernel")
  # )

  RKprime2 <- integrate(
    Vectorize(function(u) kernel_function_density_deriv2(u)^2),
    -4, 4
  )$value

  # sigma_K_2 <- switch(kernel,
  #   gaussian = 1,
  #   epanechnikov = 1 / 5,
  #   rectangular = 1 / 3,
  #   triangular = 1 / 6,
  #   biweight = 1 / 7,
  #   cosine = (pi^2 - 6) / (3 * pi^2),
  #   optcosine = (pi^2 - 8) / (pi^2),
  #   stop("unknown kernel")
  # )

  sigma_K_2 <- integrate(
    Vectorize(function(u) u^2 * kernel_function_density(u)),
    -4, 4
  )$value

  # intudW2 <- switch(kernel,
  #   gaussian = 1 / sqrt(pi),
  #   epanechnikov = 9 / 35,
  #   rectangular = 1 / 3,
  #   triangular = 7 / 30,
  #   biweight = 50 / 231,
  #   cosine = (4 * pi^2 - 15) / (12 * pi^2),
  #   optcosine = 1 / 4,
  #   stop("unknown kernel")
  # )

  intudW2 <- integrate(
    Vectorize(function(u) 2 * u * kernel_function_density(u) * kernel_function_distribution(u)),
    -4, 4
  )$value

  list(
    kernel_function_density = kernel_function_density,
    kernel_function_distribution = kernel_function_distribution,
    kernel_function_density_deriv = kernel_function_density_deriv,
    kernel_function_density_deriv2 = kernel_function_density_deriv2,
    kernel_function_conv = kernel_function_conv,
    RK = RK,
    RKprime = RKprime,
    RKprime2 = RKprime2,
    sigma_K_2 = sigma_K_2,
    intudW2 = intudW2
  )
}


.simpsons_rule <- function(x, fx) {
  # https://scicomp.stackexchange.com/questions/25649/composite-simpsons-rule-with-odd-intervals
  valid <- which(!is.na(x) & !is.na(fx) & is.finite(x) & is.finite(fx))
  x <- x[valid]
  fx <- fx[valid]

  ord <- order(x)
  x <- x[ord]
  fx <- fx[ord]
  n <- length(x)
  if (n < 5) stop("At least 5 points are required for Simpson's rule")

  h <- (x[2] - x[1])
  if (any(abs(diff(x) - h) > .Machine$double.eps^0.5)) stop("x must be equally spaced")

  integral <- fx[1] + fx[n] + 4 * sum(fx[seq.int(2, n - 1, by = 2)]) + 2 * sum(fx[seq.int(3, n - 2, by = 2)])
  integral * h / 3
}


.validate_number <- function(x, name) {
  if (!is.numeric(x) || any(is.na(x)) || any(!is.finite(x)) || length(x) != 1) {
    stop(paste("Invalid", name, ":", name, "must be numeric of lenght one"))
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

.get_prob_grid <- function(x, from, to, nb, bw, adjust) {
  if (missing(x)) {
    .validate_number(nb, "nb")
    .validate_number(from, "from")
    .validate_number(to, "to")
    if (from < 0 || to > 1) stop("'from' and 'to' must be numeric on [0,1]")
    if (from > to) stop("'from' must be smaller than 'to'")
    if (!all.equal(nb, as.integer(nb)) || nb <= 0) stop("'nb' must be a positive integer")
    x <- seq.int(from, to, length.out = nb)
  } else {
    if (!is.numeric(x) || !is.vector(x) || !all(x >= 0 & x <= 1)) stop("'x' must be a numeric vector on [0,1]")
    x <- sort(x)
  }

  message(sprintf("Interval for Estimation: [%f, %f]", min(x), max(x)))

  .validate_number(bw, "bw")
  .validate_number(adjust, "adjust")

  bw <- adjust * bw
  list(x = x, from = from, to = to, bw = bw)
}

.get_xaxn_grid <- function(y, x, from, to, nb, plot) {
  if (missing(x)) {
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
    x <- seq.int(from, to, length.out = nb)
  } else {
    if (!is.numeric(x) || !is.vector(x)) stop("'x' must be a numeric vector")
    from <- min(x)
    to <- max(x)
  }

  message(sprintf("Interval for Estimation: [%f, %f]", min(x), max(x)))

  if (!is.logical(plot)) {
    stop("argument 'plot' must be logical")
  }
  list(x = x, from = from, to = to)
}
