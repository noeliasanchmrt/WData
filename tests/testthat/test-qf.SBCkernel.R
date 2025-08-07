test_that("qf.SBCkernel() returns a valid estimate", {
  result <- qf.SBCkernel(biased_models[[1]], bw = 0.25, plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$x, "double")
  expect_type(result$est_values, "double")
  expect_gt(length(result$x), 0)
  expect_true(all(!is.na(result$est_values)))
})

test_that("qf.SBCkernel() correctly handles different bandwidth selection methods", {
  bw_methods <- list("bw.q.SBCcv", 0.25)

  results <- lapply(bw_methods, function(bw_method) {
    suppressWarnings(qf.SBCkernel(biased_models[[1]], bw = bw_method, plot = FALSE))
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


test_that("qf.SBCkernel() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    qf.SBCkernel(biased_models[[1]], kernel = k, bw = 0.25, plot = FALSE)
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


test_that("qf.SBCkernel() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("qf.SBCkernel_model_", i, "_kernel_", k),
        function() {
          qf.SBCkernel(biased_models[[i]], bw = 0.25, kernel = k, plot = TRUE)
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
})
