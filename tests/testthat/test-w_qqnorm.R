test_that("w_qqnorm() returns expected structure", {
  lapply(seq_along(biased_models), function(i) {
    result <- w_qqnorm(biased_models[[i]], plot.it = FALSE)

    expect_type(result, "list")
    expect_named(result, c("x", "y"))
    expect_type(result$x, "double")
    expect_type(result$y, "double")
    expect_length(result$x, length(biased_models[[i]]))
    expect_length(result$y, length(biased_models[[i]]))
  })
})

test_that("w_qqnorm() and w_qqline() produce stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("w_qqnorm_model_", i),
      function() {
        w_qqnorm(biased_models[[i]], plot.it = TRUE)
        w_qqline(biased_models[[i]], qmethod = "qf.SBC", col = "magenta", lty = 1)
        w_qqline(biased_models[[i]], qmethod = "qf.sen", col = "magenta", lty = 2)
      }
    )
  })
})
