test_that("qf.SBC() returns a valid estimate", {
  cdf <- sp.akbaripi(biased_models[[1]])

  expect_s3_class(cdf, "eqf")
  expect_s3_class(cdf, "stepfun")
  expect_type(cdf(0.5), "double")
})

test_that("sp.akbaripi() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("sp.akbaripi_model_", i),
      function() {
        plot(sp.akbaripi(biased_models[[i]]))
        suppressWarnings(curve(
          {
            \(.)  sp_list[[i]](.)
          }(x),
          from = 0.025,
          to = 0.975,
          col = "magenta",
          add = TRUE
        ))
      }
    )
  })
})
