test_that("qf.sen() returns a valid estimate", {
  cdf <- qf.sen(biased_models[[1]])

  expect_s3_class(cdf, "eqf")
  expect_s3_class(cdf, "stepfun")
  expect_type(cdf(0.5), "double")
  expect_gte(cdf(0), min(biased_models[[1]]))
  expect_lte(cdf(1), max(biased_models[[1]]))
})

test_that("qf.sen() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("qf.sen_model_", i),
      function() {
        plot(qf.sen(biased_models[[i]]),
          xlab = "", ylab = "",
          main = paste("Modelo", i)
        )
        suppressWarnings(curve(
          {
            \(.)  qf_list[[i]](.)
          }(x),
          col = "magenta",
          add = TRUE
        ))
      }
    )
  })
})
