test_that("bw.F.BD() returns valid bandwidth values", {
  x_vals <- seq(0.01, 1, length.out = 512)
  c_vals <- 0.25

  bw <- bw.F.BD(biased_models[[1]], y.seq = x_vals, cy.seq = c_vals)

  expect_type(bw, "double")
  expect_length(bw, length(x_vals))
  expect_true(all(bw > 0, na.rm = TRUE))
})

test_that("bw.F.BD() produces different values for different cy.seq", {
  x_vals <- seq(0.01, 1, length.out = 512)
  bw1 <- bw.F.BD(biased_models[[1]], y.seq = x_vals, cy.seq = 0.25)
  bw2 <- bw.F.BD(biased_models[[1]], y.seq = x_vals, cy.seq = 1.3)

  expect_false(identical(bw1, bw2), info = "Different cy.seq values should produce different bandwidths")
})

test_that("bw.F.BD() handles incorrect inputs correctly", {
  x_vals <- seq(0.01, 1, length.out = 512)
  c_vals <- rep(0.25, 512)

  expect_error(bw.F.BD(numeric(0), y.seq = x_vals, cy.seq = c_vals), "need at least 2 data points")
  expect_error(bw.F.BD(rep(NA, 10), y.seq = x_vals, cy.seq = c_vals), "invalid 'y'")
  expect_error(bw.F.BD(biased_models[[1]], y.seq = NULL, cy.seq = c_vals), "argument 'y.seq' must be a vector")
  expect_error(bw.F.BD(biased_models[[1]], y.seq = x_vals, cy.seq = rep(-1, 512)), "argument 'cy.seq' must be positive")
})


test_that("cdf.bd() with bw.F.BD() produces stable plots for all kernels", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.BD_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]], plot = TRUE, bw = "bw.F.BD", kernel = k, cy.seq = 0.25)
          suppressWarnings(curve(
            {
              \(.)  cdf_list[[i]](.)
            }(x),
            col = "magenta",
            add = TRUE
          ))
        }
      )
    })
  })
})
