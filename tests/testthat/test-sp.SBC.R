test_that("sp.SBC() returns a valid estimate", {
  result <- sp.SBC(biased_models[[1]], bw = 0.1, plot = FALSE)

  # expect_s3_class(result, "density")
  expect_type(result$x, "double")
  expect_type(result$est_values, "double")
  expect_gt(length(result$x), 0)
  expect_true(all(!is.na(result$est_values)))
})

test_that("sp.SBC() correctly handles different kernels", {
  results <- lapply(kernels, function(k) {
    sp.SBC(biased_models[[1]], kernel = k, bw = 0.1, plot = FALSE)
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


test_that("sp.SBC() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("sp.SBC_model_", i, "_kernel_", k),
        function() {
          sp.SBC(biased_models[[i]], bw = 0.1, kernel = k, plot = TRUE)
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
})
