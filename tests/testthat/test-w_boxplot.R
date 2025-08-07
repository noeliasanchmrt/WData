test_that("w_boxplot() handles valid inputs", {
  result <- w_boxplot(biased_models[[1]], plot = FALSE, qmethod = "qf.SBC")

  expect_type(result, "list")
  expect_named(result, c("stats", "n", "conf", "out", "group", "names"))

  expect_type(result$stats, "double")
  expect_equal(dim(result$stats), c(5, 1))

  expect_type(result$n, "double")
  expect_gt(result$n, 0)

  expect_type(result$conf, "double")
  expect_equal(dim(result$conf), c(2, 1))

  expect_type(result$out, "double")
})

test_that("w_boxplot.formula() handles valid inputs", {
  df <- data.frame(
    Width = biased_models[[1]],
    Group = sample(c("A", "B"), length(biased_models[[1]]), replace = TRUE)
  )

  result <- w_boxplot(Width ~ Group, data = df, plot = FALSE, qmethod = "qf.SBC")
  expect_type(result, "list")
  expect_named(result, c("stats", "n", "conf", "out", "group", "names"))
})

test_that("w_boxplot() produces stable plots", {
  skip_on_os(os = c("windows", "linux"))

  lapply(seq_along(biased_models), function(i) {
    vdiffr::expect_doppelganger(
      paste0("w_boxplot_qf_SBC_model_", i),
      function() {
        w_boxplot(biased_models[[i]],
          qmethod = "qf.SBC", horizontal = FALSE,
          col = rgb(0, 0, 1, 0.25), width = 0.35,
          ylab = "", pch = 16, border = c("blue")
        )

        suppressWarnings(abline(h = qf_list[[i]](c(0.25, 0.5, 0.75)), col = "magenta"))
      }
    )

    vdiffr::expect_doppelganger(
      paste0("w_boxplot_qf_sen_model_", i),
      function() {
        w_boxplot(biased_models[[i]],
          qmethod = "qf.sen", horizontal = FALSE,
          col = rgb(0, 0, 1, 0.25), width = 0.35,
          ylab = "", pch = 16, border = c("blue")
        )

        suppressWarnings(abline(h = qf_list[[i]](c(0.25, 0.5, 0.75)), col = "magenta"))
      }
    )
  })
})
