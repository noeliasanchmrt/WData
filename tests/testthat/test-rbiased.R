test_that("rbiased() generates valid samples", {
  set.seed(123) # Fix seed for reproducibility
  samples <- lapply(models_params, function(p) {
    suppressWarnings(rbiased(
      n = 500, fx = "mixnorm",
      mean = p$mean, sd = p$sd, pro = p$pro, lim = 0.01, plot = FALSE, stop = FALSE
    ))
  })

  lapply(samples, function(data) {
    expect_type(data, "double")
    expect_true(all(!is.na(data)))
    expect_gt(length(data), 0)
  })
})

test_that("rbiased() errors when stop = TRUE", {
  set.seed(123) # Fix seed for reproducibility
  lapply(models_params, function(p) {
    expect_error(
      suppressWarnings(
        rbiased(
          n = 500, fx = "mixnorm",
          mean = p$mean - shift, sd = p$sd, pro = p$pro, lim = 0.01, plot = FALSE, stop = TRUE
        ),
        regexp = "bias function"
      )
    )
  })
})

test_that("rbiased() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  set.seed(123) # Fix seed for reproducibility
  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      title = paste0("rbiased_model_", i),
      fig = function() {
        suppressWarnings(rbiased(
          n = 500, w = function(y) ifelse(y >= 0, y, NA), fx = "mixnorm",
          mean = models_params[[i]]$mean,
          sd = models_params[[i]]$sd,
          pro = models_params[[i]]$pro,
          lim = 0.01, plot = TRUE, stop = FALSE
        ))
      }
    )
  })
})
