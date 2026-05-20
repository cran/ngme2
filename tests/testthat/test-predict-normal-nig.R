test_that("predict handles normal-nig latent noise", {
  withr::local_seed(1)

  n <- 12
  data <- data.frame(
    y = as.numeric(arima.sim(list(ar = 0.4), n = n)),
    idx = seq_len(n)
  )

  fit <- ngme(
    y ~ 0 + f(idx, model = ar1(), noise = noise_normal_nig(), name = "latent"),
    data = data,
    control_opt = control_opt(
      iterations = 30,
      n_batch = 3,
      n_min_batch = 1,
      burnin = 1,
      verbose = FALSE,
      print_check_info = FALSE,
      store_traj = FALSE,
      seed = 1
    )
  )

  pred_latent <- predict(
    fit,
    map = list(latent = data$idx),
    type = "latent",
    estimator = c("mean", "sd", "0.025q", "0.975q"),
    sampling_size = 10,
    burnin_size = 2,
    seed = 1
  )
  pred_response <- predict(
    fit,
    map = list(latent = data$idx),
    type = "response",
    estimator = c("mean", "sd"),
    sampling_size = 10,
    burnin_size = 2,
    seed = 1
  )
  pred_latent_samples <- predict(
    fit,
    map = list(latent = data$idx),
    type = "latent",
    estimator = c("mean", "sd", "0.025q", "0.975q"),
    sampling_size = 10,
    burnin_size = 2,
    seed = 1,
    return_samples = TRUE
  )

  expect_length(pred_latent$mean, n)
  expect_length(pred_response$mean, n)
  expect_length(pred_response$sd, n)
  expect_true(all(is.finite(pred_latent$mean)))
  expect_true(all(is.finite(pred_response$mean)))
  expect_true(all(is.finite(pred_response$sd)))
  expect_equal(pred_latent$mean, pred_latent_samples$mean)
  expect_equal(pred_latent$sd, pred_latent_samples$sd)
  expect_equal(pred_latent$`0.025q`, pred_latent_samples$`0.025q`)
  expect_equal(pred_latent$`0.975q`, pred_latent_samples$`0.975q`)
})
