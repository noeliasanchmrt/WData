test_that("bw.q.SBCcv() returns a valid bandwidth", {
  bw <- bw.q.SBCcv(biased_models[[1]], nh = 25L, plot = FALSE)

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.q.SBCcv() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.q.SBCcv(biased_models[[1]],
      nh = 5L, kernel = k, plot = FALSE
    )
  })

  expect_false(all(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("bw.q.SBCcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw.q.SBCcv_model_", i, "_kernel_", k),
        function() {
          bw.q.SBCcv(biased_models[[i]],
            nh = 25L, kernel = k, plot = TRUE
          )
        }
      )
    })
  })
})

test_that("qf.SBCkernel() with bw.q.SBCcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("qf.SBCkernel_bw.f.BGMcv_model_", i, "_kernel_", k),
        function() {
          qf.SBCkernel(biased_models[[i]],
            bw = "bw.q.SBCcv",
            nh = 25L,
            kernel = k, plot = TRUE
          )
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
