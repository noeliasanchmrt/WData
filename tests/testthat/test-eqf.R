test_that("eqf() creates a valid eqf object", {
  obj <- qf.sen(biased_models[[1]])

  expect_s3_class(obj, "eqf")
})

test_that("print.eqf() outputs the expected format", {
  obj <- qf.sen(biased_models[[1]])
  expect_snapshot_output(print(obj))
})

test_that("summary.eqf() returns a valid summary", {
  obj <- qf.sen(biased_models[[1]])
  summary_obj <- summary(obj)

  expect_s3_class(summary_obj, "summary.ecdf")
  expect_named(summary_obj, c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
})

test_that("quantile.eqf() computes correct quantiles", {
  obj <- qf.sen(biased_models[[1]])
  q_vals <- quantile(obj, probs = c(0.25, 0.5, 0.75))

  expect_type(q_vals, "double")
  expect_length(q_vals, 3)
})

test_that("plot.eqf() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  obj <- qf.sen(biased_models[[1]])

  vdiffr::expect_doppelganger(
    "plot_eqf",
    function() plot(obj, main = "Empirical QF")
  )
})
