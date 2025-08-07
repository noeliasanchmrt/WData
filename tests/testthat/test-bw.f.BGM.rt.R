test_that("bw.f.BGM.rt() returns valid bandwidths", {
  bw <- bw.f.BGM.rt(biased_models[[1]])

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.f.BGM.rt() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.f.BGM.rt(biased_models[[1]], kernel = k)
  })

  expect_false(any(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("df.jones() with bw.f.BGM.rt() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGM.rt_model_", i, "_kernel_", k),
        function() {
          df.jones(biased_models[[i]], bw = "bw.f.BGM.rt", kernel = k, plot = TRUE)
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
