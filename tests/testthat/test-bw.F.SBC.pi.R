test_that("bw.F.SBC.pi() correctly handles pilot bandwidth selection and returns a valid bandwidth", {
  bw <- bw.F.SBC.pi(biased_models[[1]])
  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.F.SBC.pi() correctly handles different kernels", {
  supported_kernels <- c("gaussian", "epanechnikov", "biweight", "cosine", "optcosine")

  bw_values <- sapply(supported_kernels, function(k) {
    bw.F.SBC.pi(biased_models[[1]], kernel = k)
  })

  expect_false(any(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("bw.F.SBC.pi() throws an error for unsupported kernels", {
  unsupported_kernels <- c("rectangular")

  lapply(unsupported_kernels, function(k) {
    expect_error(
      bw.F.SBC.pi(biased_models[[1]], kernel = k),
      "not supported for plug-in pilot bandwidth"
    )
  })
})

test_that("cdf.bd() with bw.F.SBC.pi() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  supported_kernels <- c("gaussian", "biweight", "cosine", "optcosine")

  lapply(seq_along(biased_models), function(i) {
    lapply(supported_kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("cdf.bd_bw.F.SBC.pi_model_", i, "_kernel_", k),
        function() {
          cdf.bd(biased_models[[i]], bw = "bw.F.SBC.pi", kernel = k, plot = TRUE)
          suppressWarnings(curve(
            {
              \(.) cdf_list[[i]](.)
            }(x),
            col = "magenta",
            add = TRUE
          ))
        }
      )
    })
  })
})
