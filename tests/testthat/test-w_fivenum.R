test_that("w_fivenum() returns a valid five-number summary", {
  lapply(biased_models, function(data) {
    summary <- w_fivenum(data)

    expect_type(summary, "double")
    expect_length(summary, 5)
    expect_true(all(!is.na(summary)))
    expect_lte(summary[1], summary[2])
    expect_lte(summary[2], summary[3])
    expect_lte(summary[3], summary[4])
    expect_lte(summary[4], summary[5])
  })
})

test_that("w_fivenum() correctly handles different quantile methods", {
  lapply(biased_models, function(data) {
    summary_SBC <- w_fivenum(data, qmethod = "qf.SBC")
    summary_sen <- w_fivenum(data, qmethod = "qf.sen")

    expect_length(summary_SBC, 5)
    expect_length(summary_sen, 5)
    expect_false(identical(summary_SBC, summary_sen),
      info = "Different quantile methods should produce different results"
    )
  })
})

test_that("w_fivenum() correctly handles missing values", {
  data_with_na <- c(biased_models[[1]], NA, NA)

  summary_rm <- w_fivenum(data_with_na, na.rm = TRUE)
  summary_narm <- w_fivenum(data_with_na, na.rm = FALSE)

  expect_length(summary_rm, 5)
  expect_true(all(!is.na(summary_rm)))

  expect_length(summary_narm, 5)
  expect_true(all(is.na(summary_narm)))
})

test_that("w_fivenum() throws an error for invalid inputs", {
  expect_error(
    w_fivenum(biased_models[[1]], qmethod = "invalid_method"),
    "'arg' should be"
  )
})
