test_that("df.jones() returns a valid estimate", {
  result <- df.jones(biased_models[[1]], bw = "bw.f.BGMnrd0", plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$x, "double")
  expect_type(result$est_values, "double")
  expect_gt(length(result$x), 0)
  expect_true(all(!is.na(result$est_values)))
})

test_that("df.jones() correctly handles different bandwidth selection methods", {
  bw_methods <- list("bw.f.BGMnrd0", "bw.f.BGMcv", "bw.f.BGMboot1", "bw.f.BGMboot2", 0.5)

  results <- lapply(bw_methods, function(bw_method) {
    suppressWarnings(df.jones(biased_models[[1]], bw = bw_method, plot = FALSE))
  })

  for (i in seq_along(results)) {
    for (j in seq_along(results)) {
      if (i != j) {
        expect_false(identical(results[[i]]$est_values, results[[j]]$est_values),
          info = paste("Bandwidth comparison failed for", bw_methods[i], "vs", bw_methods[j])
        )
      }
    }
  }
})

test_that("df.jones() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    df.jones(biased_models[[1]], kernel = k, bw = "bw.f.BGMnrd0", plot = FALSE)
  })

  for (i in seq_along(results)) {
    for (j in seq_along(results)) {
      if (i != j) {
        expect_false(identical(results[[i]]$est_values, results[[j]]$est_values),
          info = paste("Kernel comparison failed for", kernels[i], "vs", kernels[j])
        )
      }
    }
  }
})
