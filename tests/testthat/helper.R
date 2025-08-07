options(vdiffr.graphics.engine = "svglite")

kernels <- c(
  "gaussian",
  "epanechnikov",
  "rectangular",
  "triangular",
  "biweight",
  "cosine",
  "optcosine"
)

shift <- 4

models_params <- list(
  m1 = list(mean = 0 + shift, sd = 1, pro = 1),
  m2 = list(mean = c(0, 1 / 2, 12 / 13) + shift, sd = c(1, 2 / 3, 5 / 9), pro = c(1 / 5, 1 / 5, 3 / 5)),
  m3 = list(mean = 3 * ((2 / 3)^(0:7) - 1) + shift, sd = (2 / 3)^(0:7), pro = rep(1 / 8, 8)),
  m4 = list(mean = c(0, 0) + shift, sd = c(1, 1 / 10), pro = c(2 / 3, 1 / 3)),
  m5 = list(mean = c(-3 / 2, 3 / 2) + shift, sd = c(1 / 2, 1 / 2), pro = c(1 / 2, 1 / 2)),
  m6 = list(mean = c(0, 3 / 2) + shift, sd = c(1, 1 / 3), pro = c(3 / 4, 1 / 4)),
  m7 = list(mean = c(0, (0:4) / 2 - 1) + shift, sd = c(1, rep(1 / 10, 5)), pro = c(1 / 2, rep(1 / 10, 5))),
  m8 = list(mean = (65 - 96 * (1 / 2)^(0:5)) / 21 + shift, sd = (32 / 63) * (1 / 2)^(0:5), pro = (2^(5 - (0:5))) / 63)
)

create_density <- function(mean, sd, pro) {
  function(x) KScorrect::dmixnorm(x, mean = mean, sd = sd, pro = pro)
}

create_cdf <- function(mean, sd, pro) {
  function(x) KScorrect::pmixnorm(x, mean = mean, sd = sd, pro = pro)
}

create_quantile <- function(mean, sd, pro) {
  function(x) KScorrect::qmixnorm(x, mean = mean, sd = sd, pro = pro)
}

create_derivative <- function(mean, sd, pro) {
  function(x) sapply(x, function(xi) sum(-pro * dnorm(xi, mean = mean, sd = sd) * (xi - mean) / (sd^2)))
}

df_list <- lapply(models_params, function(p) create_density(p$mean, p$sd, p$pro))
cdf_list <- lapply(models_params, function(p) create_cdf(p$mean, p$sd, p$pro))
qf_list <- lapply(models_params, function(p) create_quantile(p$mean, p$sd, p$pro))
sp_list <- Map(function(f, g) function(x) 1 / f(g(x)), df_list, qf_list)

set.seed(1234)
biased_models <- lapply(models_params, function(p) {
  rbiased(
    n = 100, w = function(y) ifelse(y >= 0, y, NA), fx = "mixnorm",
    mean = p$mean, sd = p$sd, pro = p$pro, lim = 0.01, plot = FALSE, stop = FALSE
  )
})

biased_models <- biased_models[1]
