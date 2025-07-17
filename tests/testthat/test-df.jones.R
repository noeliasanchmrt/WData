test_that("df.jones() returns a valid estimate", {
  result <- df.jones(biased_models[[1]], bw = "bw.f.BGM.rt", plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$y.seq, "double")
  expect_type(result$f.hat, "double")
  expect_gt(length(result$y.seq), 0)
  expect_true(all(!is.na(result$f.hat)))
})

test_that("df.jones() correctly handles different bandwidth selection methods", {
  bw_methods <- list("bw.f.BGM.rt", "bw.f.BGM.cv", "bw.f.BGM.boot1", "bw.f.BGM.boot2", 0.5)

  results <- lapply(bw_methods, function(bw_method) {
    suppressWarnings(df.jones(biased_models[[1]], bw = bw_method, plot = FALSE))
  })

  for (i in seq_along(results)) {
    for (j in seq_along(results)) {
      if (i != j) {
        expect_false(identical(results[[i]]$f.hat, results[[j]]$f.hat),
          info = paste("Bandwidth comparison failed for", bw_methods[i], "vs", bw_methods[j])
        )
      }
    }
  }
})

test_that("df.jones() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    df.jones(biased_models[[1]], kernel = k, bw = "bw.f.BGM.rt", plot = FALSE)
  })

  for (i in seq_along(results)) {
    for (j in seq_along(results)) {
      if (i != j) {
        expect_false(identical(results[[i]]$f.hat, results[[j]]$f.hat),
          info = paste("Kernel comparison failed for", kernels[i], "vs", kernels[j])
        )
      }
    }
  }
})
