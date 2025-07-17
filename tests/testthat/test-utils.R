test_that(".check_biased_sample() handles valid inputs", {
  result <- .check_biased_sample(1:10, function(y) y + 1)

  expect_type(result, "list")
  expect_named(result, c("data.name", "n", "y", "w", "yw", "weights", "uw", "uwb", "has.na"))
  expect_equal(result$n, 10)
  expect_type(result$w, "closure")
  expect_length(result$y, 10)
  expect_length(result$yw, 10)
})

test_that(".check_biased_sample() throws errors for invalid inputs", {
  expect_error(.check_biased_sample("not_numeric", function(y) y + 1), "invalid 'y'")
  expect_error(.check_biased_sample(c(1, NA), function(y) y + 1), "need at least 2 data points")
  expect_error(.check_biased_sample(1:10, "not_a_function"), "argument 'w' must be a function")
  expect_error(.check_biased_sample(1:10, function(y) NA), "function 'w' must be bounded")
})

test_that(".get_kernel_values() handles valid inputs", {
  kernels <- c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")

  lapply(kernels, function(k) {
    result <- .get_kernel_values(k)

    expect_type(result, "list")
    expect_named(result, c(
      "kernel_function_density",
      "kernel_function_distribution",
      "kernel_function_density_deriv1",
      "kernel_function_density_deriv2",
      "kernel_function_conv",
      "kernel_r",
      "kernel_r_deriv1",
      "kernel_r_deriv2",
      "kernel_eta",
      "kernel_kappa"
    ))
    expect_type(result$kernel_function_density, "closure")
  })
})

test_that(".get_kernel_values() throws errors for invalid inputs", {
  expect_error(.get_kernel_values("unknown_kernel"), "unknown kernel")
})

test_that(".get_kernel_values() is consistent", {
  kernels <- c(
    "gaussian", "epanechnikov", "rectangular",
    "triangular",
    "biweight", "cosine", "optcosine"
  )

  for (i in kernels) {
    k <- .get_kernel_values(i)

    expect_lt(abs(integrate(
      Vectorize(function(u) k$kernel_function_density(u)^2),
      -4, 4
    )$value - k$kernel_r), 0.01)

    expect_lt(abs(integrate(
      Vectorize(function(u) k$kernel_function_density_deriv1(u)^2),
      -4, 4
    )$value - k$kernel_r_deriv1), 0.01)

    expect_lt(abs(integrate(
      Vectorize(function(u) k$kernel_function_density_deriv2(u)^2),
      -4, 4
    )$value - k$kernel_r_deriv2), 0.01)

    expect_lt(abs(integrate(
      Vectorize(function(u) u^2 * k$kernel_function_density(u)),
      -4, 4
    )$value - k$kernel_eta), 0.01)

    expect_lt(abs(integrate(
      Vectorize(function(u) 2 * u * k$kernel_function_density(u) * k$kernel_function_distribution(u)),
      -4, 4
    )$value - k$kernel_kappa), 0.01)
  }
})


test_that(".simpsons_rule() handles valid inputs", {
  x <- seq(0, 1, length.out = 11)
  fx <- x^2

  result <- .simpsons_rule(x, fx)
  expect_type(result, "double")
  expect_true(abs(result - 1 / 3) < 1e-3) # Expected integral of x^2 from 0 to 1 is 1/3
})

test_that(".validate_number() handles valid inputs", {
  expect_silent(.validate_number(1, "test"))
  expect_silent(.validate_number(3.14, "test"))
  expect_silent(.validate_number(-5, "test"))
})

test_that(".validate_number() throws errors for invalid inputs", {
  expect_error(.validate_number(NA, "test"), "Invalid test")
  expect_error(.validate_number(NaN, "test"), "Invalid test")
  expect_error(.validate_number(Inf, "test"), "Invalid test")
  expect_error(.validate_number("string", "test"), "Invalid test")
  expect_error(.validate_number(c(1, 2), "test"), "Invalid test")
  expect_error(.validate_number(list(1), "test"), "Invalid test")
})

test_that(".get_bandwidth_grid() handles valid inputs", {
  grid <- .get_bandwidth_grid(nh = 5, lower = 0.1, upper = 1, tol = 0.01, plot = FALSE)

  expect_type(grid, "double")
  expect_length(grid, 5)
  expect_true(all(grid >= 0.1 & grid <= 1))
})

test_that(".get_bandwidth_grid() throws errors for invalid inputs", {
  expect_error(.get_bandwidth_grid(nh = 0, lower = 0.1, upper = 1, tol = 0.01, plot = FALSE), "'nh' must be a positive integer")
  expect_error(.get_bandwidth_grid(nh = 5, lower = -1, upper = 1, tol = 0.01, plot = FALSE), "'from' and 'to' must be positive numbers")
  expect_error(.get_bandwidth_grid(nh = 5, lower = 0.5, upper = 0.1, tol = 0.01, plot = FALSE), "'lower' must be smaller than 'upper'")
})

test_that(".get_xaxn_grid() handles valid inputs", {
  result <- .get_xaxn_grid(1:10, from = 1, to = 10, nb = 5, plot = FALSE)

  expect_type(result, "list")
  expect_named(result, c("y.seq", "from", "to"))
  expect_length(result$y.seq, 5)
  expect_equal(result$from, 1)
  expect_equal(result$to, 10)
})

test_that(".get_xaxn_grid() throws errors for invalid inputs", {
  expect_error(.get_xaxn_grid(1:10, y.seq = NULL, from = 10, to = 1, nb = 5, plot = FALSE), "'y.seq' must be a numeric vector")
  expect_error(.get_xaxn_grid(1:10, from = 10, to = 1, nb = 5, plot = FALSE), "'from' must be smaller than 'to'")
  expect_error(.get_xaxn_grid(1:10, from = 1, to = 10, nb = -5, plot = FALSE), "'nb' must be a positive integer")
})
