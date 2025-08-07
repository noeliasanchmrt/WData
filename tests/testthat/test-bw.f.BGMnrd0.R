test_that("bw.f.BGMnrd0() returns valid bandwidths", {
  bw <- bw.f.BGMnrd0(biased_models[[1]])

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.f.BGMnrd0() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.f.BGMnrd0(biased_models[[1]], kernel = k)
  })

  expect_false(any(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("df.jones() with bw.f.BGMnrd0() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGMnrd0_model_", i, "_kernel_", k),
        function() {
          df.jones(biased_models[[i]], bw = "bw.f.BGMnrd0", kernel = k, plot = TRUE)
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
