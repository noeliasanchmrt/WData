test_that("df.bhatta() returns a valid estimate", {
  result <- df.bhatta(biased_models[[1]], plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$y.seq, "double")
  expect_type(result$f.hat, "double")
  expect_true(all(!is.na(result$f.hat)))
  expect_gt(length(result$y.seq), 0)
})

test_that("df.bhatta() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    df.bhatta(biased_models[[1]], bw = "nrd0", kernel = k, plot = FALSE)
  })

  # Ensure all kernel results are different
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

test_that("df.bhatta() handles different bandwidth selection methods", {
  result1 <- df.bhatta(biased_models[[1]], bw = "nrd0", plot = FALSE)
  result2 <- df.bhatta(biased_models[[1]], bw = "ucv", plot = FALSE)

  expect_false(identical(result1$f.hat, result2$f.hat))
})

test_that("df.bhatta() produces stable plots for all kernels", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.bhatta_model_", i, "_kernel_", k),
        function() {
          df.bhatta(biased_models[[i]], kernel = k, plot = TRUE)
          suppressWarnings(curve(
            {
              \(.)  df_list[[i]](.)
            }(x),
            col = "magenta",
            add = TRUE
          ))
        }
      )
    })
  })
})
