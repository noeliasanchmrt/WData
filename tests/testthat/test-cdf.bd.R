test_that("cdf.bd() returns a valid estimate", {
  result <- cdf.bd(biased_models[[1]], bw = "bw.F.SBCnrd0", plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$x, "double")
  expect_type(result$est_values, "double")
  expect_gt(length(result$x), 0)
  expect_true(all(!is.na(result$est_values)))
})

test_that("cdf.bd() correctly handles different bandwidth selection methods", {
  bw_methods <- list("bw.F.SBCnrd0", "bw.F.SBCcv", "bw.F.SBCpi", 0.5)

  results <- lapply(bw_methods, function(bw_method) {
    cdf.bd(biased_models[[1]], bw = bw_method, plot = FALSE)
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

test_that("cdf.bd() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    cdf.bd(biased_models[[1]], kernel = k, bw = "bw.F.SBCnrd0", plot = FALSE)
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


test_that("cdf.bd() with boundary correction produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]], plot = TRUE, kernel = k)
          both <- cdf.bd(biased_models[[i]], kernel = k, correction = "both", plot = F)
          right <- cdf.bd(biased_models[[i]], kernel = k, correction = "right", plot = F)
          left <- cdf.bd(biased_models[[i]], kernel = k, correction = "left", plot = F)
          lines(both$x, both$est_values, col = "blue", lty = 2)
          lines(right$x, right$est_values, col = "blue", lty = 3)
          lines(left$x, left$est_values, col = "blue", lty = 4)
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
})
