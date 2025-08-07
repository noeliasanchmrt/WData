test_that("cdf.cox() returns a valid estimate", {
  cdf <- cdf.cox(biased_models[[1]])

  expect_s3_class(cdf, "ecdf")
  expect_s3_class(cdf, "stepfun")
  expect_type(cdf(0.5), "double")
  expect_gte(cdf(min(biased_models[[1]])), 0)
  expect_lte(cdf(max(biased_models[[1]])), 1)
})

test_that("cdf.cox() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("cdf.cox_model_", i),
      function() {
        plot(cdf.cox(biased_models[[i]]))
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
