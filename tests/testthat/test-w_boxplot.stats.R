test_that("w_boxstats() returns valid boxplot statistics", {
  lapply(biased_models, function(data) {
    result <- w_boxstats(data)

    expect_type(result$stats, "double")
    expect_length(result$stats, 5)
    expect_gte(result$n, 1)
    expect_lte(result$stats[1], result$stats[2])
    expect_lte(result$stats[2], result$stats[3])
    expect_lte(result$stats[3], result$stats[4])
    expect_lte(result$stats[4], result$stats[5])
  })
})

test_that("w_boxstats() correctly handles different quantile methods", {
  lapply(biased_models, function(data) {
    result_SBC <- w_boxstats(data, qmethod = "qf.SBC")
    result_sen <- w_boxstats(data, qmethod = "qf.sen")

    expect_length(result_SBC$stats, 5)
    expect_length(result_sen$stats, 5)
    expect_false(identical(result_SBC$stats, result_sen$stats),
      info = "Different quantile methods should produce different results"
    )
  })
})

test_that("w_boxstats() correctly handles coef = 0", {
  result <- w_boxstats(biased_models[[1]], coef = 0)

  expect_equal(result$stats[1], min(biased_models[[1]], na.rm = TRUE))
  expect_equal(result$stats[5], max(biased_models[[1]], na.rm = TRUE))
  expect_equal(length(result$out), 0)
})

test_that("w_boxstats() correctly handles confidence intervals", {
  result_with_conf <- w_boxstats(biased_models[[1]], do.conf = TRUE)
  result_without_conf <- w_boxstats(biased_models[[1]], do.conf = FALSE)

  expect_true(!is.null(result_with_conf$conf))
  expect_true(is.null(result_without_conf$conf))
})

test_that("w_boxstats() correctly handles outliers", {
  result_with_outliers <- w_boxstats(biased_models[[1]], do.out = TRUE)
  result_without_outliers <- w_boxstats(biased_models[[1]], do.out = FALSE)

  expect_true(is.numeric(result_with_outliers$out))
  expect_equal(length(result_without_outliers$out), 0)
})

test_that("w_boxstats() throws an error for invalid inputs", {
  expect_error(w_boxstats(biased_models[[1]], coef = -1), "'coef' must not be negative")
})
