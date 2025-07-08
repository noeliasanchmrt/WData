test_that("bw.F.SBCrt() returns valid bandwidths", {
  bw <- bw.F.SBCrt(biased_models[[1]])

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.F.SBCrt() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.F.SBCrt(biased_models[[1]], kernel = k)
  })

  expect_false(any(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("cdf.bd() with bw.F.SBCrt() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.SBCrt_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]], bw = "bw.F.SBCrt", kernel = k, plot = TRUE)
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
