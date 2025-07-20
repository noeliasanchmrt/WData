test_that("bw.f.BGM.boot2() returns a valid bandwidth", {
  bw <- suppressWarnings(bw.f.BGM.boot2(biased_models[[1]], nh = 10L, plot = FALSE))

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.f.BGM.boot2() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    suppressWarnings(bw.f.BGM.boot2(biased_models[[1]],
      nh = 10L,
      lower = IQR(biased_models[[1]]) * length(biased_models[[1]])^{
        -0.2
      } * 10^{
        -1
      },
      upper = IQR(biased_models[[1]]) *
        (log(length(biased_models[[1]])) / length(biased_models[[1]]))^{
          0.2
        } * 3, kernel = k, plot = FALSE
    ))
  })

  expect_false(all(duplicated(bw_values)), info = "Different kernels should produce different bandwidths")
})

test_that("bw.f.BGM.boot2() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw.f.BGM.boot2_model_", i, "_kernel_", k),
        function() {
          original_par <- par(no.readonly = TRUE) # Save current graphical parameters
          on.exit(par(original_par)) # Restore after test

          par(mfrow = c(1, 3)) # Ensure consistent plot layout

          suppressWarnings(bw.f.BGM.boot2(biased_models[[i]],
            kernel = k,
            from = max(0.02, min(biased_models[[i]]) - (sort(biased_models[[i]])[2] - min(biased_models[[i]]))),
            nh = 10L,
            lower = IQR(biased_models[[1]]) * length(biased_models[[1]])^{
              -0.2
            } * 10^{
              -1
            },
            upper = IQR(biased_models[[1]]) *
              (log(length(biased_models[[1]])) / length(biased_models[[1]]))^{
                0.2
              } * 3,
            plot = TRUE
          ))
        }
      )
    })
  })
})


test_that("df.jones() with bw.f.BGM.boot2() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGM.boot2_model_", i, "_kernel_", k),
        function() {
          suppressWarnings(df.jones(biased_models[[i]],
            bw = "bw.f.BGM.boot2",
            from = max(0.02, min(biased_models[[i]]) - (sort(biased_models[[i]])[2] - min(biased_models[[i]]))),
            nh = 10L,
            lower = IQR(biased_models[[1]]) * length(biased_models[[1]])^{
              -0.2
            } * 10^{
              -1
            },
            upper = IQR(biased_models[[1]]) *
              (log(length(biased_models[[1]])) / length(biased_models[[1]]))^{
                0.2
              } * 3,
            kernel = k,
            plot = TRUE
          ))
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
