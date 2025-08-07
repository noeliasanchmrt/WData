test_that("bw.f.BGMboot2() returns a valid bandwidth", {
  bw <- suppressWarnings(bw.f.BGMboot2(biased_models[[1]], nh = 25L, plot = FALSE))

  expect_type(bw, "double")
  expect_gt(bw, 0)
})

test_that("bw.f.BGMboot2() correctly handles different kernels", {
  bw_values <- sapply(kernels, function(k) {
    suppressWarnings(bw.f.BGMboot2(biased_models[[1]], nh = 25L, kernel = k, plot = FALSE))
  })

  expect_false(all(duplicated(bw_values)), info = "Different kernels should produce different bandwidths")
})

test_that("bw.f.BGMboot2() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("bw_f_BGMboot2_model_", i, "_kernel_", k),
        function() {
          original_par <- par(no.readonly = TRUE) # Save current graphical parameters
          on.exit(par(original_par)) # Restore after test

          par(mfrow = c(1, 3)) # Ensure consistent plot layout

          suppressWarnings(bw.f.BGMboot2(biased_models[[i]],
            kernel = k,
            from = max(0.02, min(biased_models[[i]]) - (sort(biased_models[[i]])[2] - min(biased_models[[i]]))),
            nh = 25L,
            plot = TRUE
          ))
        }
      )
    })
  })
})


test_that("df.jones() with bw.f.BGMboot2() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  old_width <- options("width")$width
  on.exit(options(width = old_width), add = TRUE)

  lapply(seq_along(biased_models), function(i) {
    lapply(kernels, function(k) {
      vdiffr::expect_doppelganger(
        paste0("df.jones_bw.f.BGMboot2_model_", i, "_kernel_", k),
        function() {
          suppressWarnings(df.jones(biased_models[[i]],
            bw = "bw.f.BGMboot2", kernel = k,
            from = max(0.02, min(biased_models[[i]]) - (sort(biased_models[[i]])[2] - min(biased_models[[i]]))),
            nh = 25L,
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
