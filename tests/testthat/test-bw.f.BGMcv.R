test_that("bw.f.BGMcv() returns a valid bandwidth", {
  bw <- bw.f.BGMcv(biased_models[[1]], nh = 25L, plot = FALSE)

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.f.BGMcv() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    bw.f.BGMcv(biased_models[[1]],
      lower = IQR(biased_models[[1]]) * length(biased_models[[1]])^{
        -0.2
      } * 1000^{
        -1
      },
      upper = IQR(biased_models[[1]]) *
        (log(length(biased_models[[1]])) / length(biased_models[[1]]))^{
          0.2
        } * 3,
      nh = 25L, kernel = k, plot = FALSE
    )
  })

  expect_false(all(duplicated(bw_values)),
    info = "Different kernels should produce different bandwidths"
  )
})

test_that("bw.f.BGMcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw.f.BGMcv_model_", i, "_kernel_", k),
        function() {
          original_par <- par(no.readonly = TRUE) # Save current graphical parameters
          on.exit(par(original_par)) # Restore them after the test

          par(mfrow = c(1, 3)) # Set layout for three plots

          bw.f.BGMcv(biased_models[[i]],
            lower = IQR(biased_models[[i]]) * length(biased_models[[i]])^{
              -0.2
            } * 1000^{
              -1
            },
            upper = IQR(biased_models[[i]]) *
              (log(length(biased_models[[i]])) / length(biased_models[[i]]))^{
                0.2
              } * 3,
            nh = 25L, kernel = k, plot = TRUE
          )
        }
      )
    })
  })
})

test_that("df.jones() with bw.f.BGMcv() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGMcv_model_", i, "_kernel_", k),
        function() {
          df.jones(biased_models[[i]],
            bw = "bw.f.BGMcv",
            lower = IQR(biased_models[[i]]) * length(biased_models[[i]])^{
              -0.2
            } * 1000^{
              -1
            },
            upper = IQR(biased_models[[i]]) *
              (log(length(biased_models[[i]])) / length(biased_models[[i]]))^{
                0.2
              } * 3,
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
