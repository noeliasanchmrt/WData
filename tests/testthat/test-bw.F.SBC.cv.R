test_that("bw.F.SBC.cv() returns a valid bandwidth", {
  bw <- bw.F.SBC.cv(biased_models[[1]], nh = 25L, plot = FALSE)

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.F.SBC.cv() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.F.SBC.cv(biased_models[[1]], nh = 25L, kernel = k, plot = FALSE)
  })

  expect_false(all(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("bw.F.SBC.cv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw.F.SBC.cv_model_", i, "_kernel_", k),
        function() {
          bw.F.SBC.cv(biased_models[[i]],
            nh = 25L, kernel = k, plot = TRUE
          )
        }
      )
    })
  })
})

test_that("cdf.bd() with bw.F.SBC.cv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.SBC.cv_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]],
            bw = "bw.F.SBC.cv",
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
