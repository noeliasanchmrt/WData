test_that("cdf.efro() returns a valid estimate", {
  result <- cdf.efro(biased_models[[1]], plot = FALSE)

  expect_type(result$fords, "double")
  expect_type(result$Fords, "double")
  expect_gt(length(result$x), 0)
  expect_true(all(!is.na(result$fords)))
  expect_true(all(!is.na(result$Fords)))
  expect_gte(min(result$Fords), 0)
  expect_lte(max(result$Fords), 1)
})


test_that("cdf.efro() handles errors correctly", {
  expect_error(cdf.efro(numeric(0)), "need at least 2 data points")
  expect_error(cdf.efro(rep(NA, 10)), "invalid 'y'")
  expect_error(cdf.efro(biased_models[[1]], Jn = -1), "Jn must be a positive number")
  expect_error(cdf.efro(biased_models[[1]], Ct = 0), "Ct must be a positive number")
  expect_error(cdf.efro(biased_models[[1]], Cjm = -5), "Cjm must be a positive number")
})

test_that("cdf.efro() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("cdf.efro_model_", i),
      function() {
        original_par <- par(no.readonly = TRUE) # Save current graphical parameters
        on.exit(par(original_par)) # Restore after test

        par(mfrow = c(2, 4)) # Ensure consistent plot layout

        cdf.efro(biased_models[[i]], plot = TRUE)
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
