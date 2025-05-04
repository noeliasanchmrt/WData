test_that("bw.F.SBCnrd0() returns valid bandwidths", {
  bw <- bw.F.SBCnrd0(biased_models[[1]])

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.F.SBCnrd0() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.F.SBCnrd0(biased_models[[1]], kernel = k)
  })

  expect_false(any(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("cdf.bd() with bw.F.SBCnrd0() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.SBCnrd0_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]], bw = "bw.F.SBCnrd0", kernel = k, plot = TRUE)
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
