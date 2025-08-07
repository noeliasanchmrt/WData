test_that("bw.F.SBCcv() returns a valid bandwidth", {
  bw <- bw.F.SBCcv(biased_models[[1]], nh = 25L, plot = FALSE)

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.F.SBCcv() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.F.SBCcv(biased_models[[1]], nh = 25L, kernel = k, plot = FALSE)
  })

  expect_false(all(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("bw.F.SBCcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw.F.SBCcv_model_", i, "_kernel_", k),
        function() {
          bw.F.SBCcv(biased_models[[i]],
            nh = 25L, kernel = k, plot = TRUE
          )
        }
      )
    })
  })
})

test_that("cdf.bd() with bw.F.SBCcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.SBCcv_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]],
            bw = "bw.F.SBCcv",
            nh = 25L,
            kernel = k, plot = TRUE
          )
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
