test_that("bw.f.BGMboot1() correctly handles pilot bandwidth selection and returns a valid bandwidth", {
  bw_rt <- bw.f.BGMboot1(biased_models[[1]], bw0 = "RT")
  bw_opt <- bw.f.BGMboot1(biased_models[[1]], bw0 = "Opt")

  expect_type(bw_opt, "double")
  expect_type(bw_rt, "double")

  expect_gt(bw_opt, 0)
  expect_gt(bw_rt, 0)
  expect_false(identical(bw_rt, bw_opt), info = "RT and Opt should produce different values")
})

test_that("bw.f.BGMboot1() correctly handles different kernels", {
  supported_kernels <- c("gaussian", "epanechnikov", "biweight", "cosine", "optcosine")

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  bw_values <- sapply(supported_kernels, function(k) {
    bw.f.BGMboot1(biased_models[[1]], kernel = k, bw0 = "RT")
  })

  expect_false(any(duplicated(bw_values)), info = "Different kernels should produce different bandwidths")
})

test_that("bw.f.BGMboot1() (Opt) throws an error for unsupported kernels", {
  unsupported_kernels <- c("rectangular", "triangular")

  lapply(unsupported_kernels, function(k) {
    expect_error(
      bw.f.BGMboot1(biased_models[[1]], kernel = k, bw0 = "Opt"),
      "not supported for automatic bandwidth selection"
    )
  })
})

test_that("bw.f.BGMboot1() (RT) throws an error for unsupported kernels", {
  unsupported_kernels <- c("rectangular", "triangular")

  lapply(unsupported_kernels, function(k) {
    expect_error(
      bw.f.BGMboot1(biased_models[[1]], kernel = k, bw0 = "RT"),
      "not supported for automatic bandwidth selection"
    )
  })
})

test_that("df.jones() with bw.f.BGMboot1() (RT) produces stable plots for supported kernels", {
  skip_on_os(os = c("windows", "linux"))

  supported_kernels <- c("gaussian", "epanechnikov", "biweight", "cosine", "optcosine")

  lapply(seq_along(biased_models), function(i) {
    lapply(supported_kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGMboot1_RT_model_", i, "_kernel_", k),
        function() {
          df.jones(biased_models[[i]], bw = "bw.f.BGMboot1", bw0 = "RT", kernel = k, plot = TRUE)
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


test_that("df.jones() with bw.f.BGMboot1() (Opt) produces stable plots for supported kernels", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  supported_kernels <- c("gaussian", "epanechnikov", "biweight", "cosine", "optcosine")

  lapply(seq_along(biased_models), function(i) {
    lapply(supported_kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGMboot1_Opt_model_", i, "_kernel_", k),
        function() {
          df.jones(biased_models[[i]], bw = "bw.f.BGMboot1", bw0 = "Opt", kernel = k, plot = TRUE)
          suppressWarnings(curve(
            {
              \(.) df_list[[i]](.)
            }(x),
            col = "magenta",
            add = TRUE
          ))
        }
      )
    })
  })
})
