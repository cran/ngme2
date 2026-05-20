test_that("stationary GAL moments match gamma mixing formulas", {
  moments <- noise_gal_moments(mu = 1, sigma = 2, nu = 3)

  expect_equal(names(moments), c(
    "skewness", "kurtosis", "excess_kurtosis",
    "variance", "sd", "central_moment3", "central_moment4"
  ))
  expect_false("mean" %in% names(moments))
  expect_equal(moments[["variance"]], 13 / 3)
  expect_equal(moments[["central_moment3"]], 38 / 9)
  expect_equal(moments[["central_moment4"]], 701 / 9)
  expect_equal(moments[["skewness"]], (38 / 9) / (13 / 3)^(3 / 2))
  expect_equal(moments[["kurtosis"]], (701 / 9) / (13 / 3)^2)
  expect_equal(moments[["excess_kurtosis"]], (701 / 9) / (13 / 3)^2 - 3)
})

test_that("stationary NIG moments match inverse-Gaussian mixing formulas", {
  moments <- noise_nig_moments(mu = 1, sigma = 2, nu = 3)

  expect_equal(names(moments), c(
    "skewness", "kurtosis", "excess_kurtosis",
    "variance", "sd", "central_moment3", "central_moment4"
  ))
  expect_false("mean" %in% names(moments))
  expect_equal(moments[["variance"]], 13 / 3, tolerance = 1e-12)
  expect_equal(moments[["central_moment3"]], 13 / 3, tolerance = 1e-12)
  expect_equal(moments[["central_moment4"]], 728 / 9, tolerance = 1e-12)
  expect_equal(moments[["skewness"]], (13 / 3) / (13 / 3)^(3 / 2), tolerance = 1e-12)
  expect_equal(moments[["kurtosis"]], (728 / 9) / (13 / 3)^2, tolerance = 1e-12)
  expect_equal(moments[["excess_kurtosis"]], (728 / 9) / (13 / 3)^2 - 3, tolerance = 1e-12)
})

test_that("stationary NIG moments match standard NIG closed-form formulas", {
  mu <- 1
  sigma <- 2
  nu <- 3

  moments <- noise_nig_moments(mu = mu, sigma = sigma, nu = nu)

  alpha <- sqrt(mu^2 / sigma^4 + nu / sigma^2)
  beta <- mu / sigma^2
  delta <- sqrt(nu) * sigma
  gamma <- sqrt(alpha^2 - beta^2)

  expect_equal(moments[["variance"]], delta * alpha^2 / gamma^3, tolerance = 1e-12)
  expect_equal(moments[["skewness"]], 3 * beta / sqrt(alpha^2 * delta * gamma), tolerance = 1e-12)
  expect_equal(
    moments[["excess_kurtosis"]],
    3 * (1 + 4 * beta^2 / alpha^2) / (delta * gamma),
    tolerance = 1e-12
  )
  expect_equal(
    moments[["kurtosis"]],
    3 + 3 * (1 + 4 * beta^2 / alpha^2) / (delta * gamma),
    tolerance = 1e-12
  )
})

test_that("stationary NIG and GAL object methods match direct moment functions", {
  nig_noise <- noise_nig(mu = 1, sigma = 2, nu = 3)
  gal_noise <- noise_gal(mu = 1, sigma = 2, nu = 3)

  expect_equal(noise_moments(nig_noise), noise_nig_moments(mu = 1, sigma = 2, nu = 3))
  expect_equal(noise_moments(gal_noise), noise_gal_moments(mu = 1, sigma = 2, nu = 3))
})

test_that("stationary moment object method accepts repeated stationary bases", {
  nig_noise <- update_noise(noise_nig(mu = 1, sigma = 2, nu = 3), n = 4)
  gal_noise <- update_noise(noise_gal(mu = 1, sigma = 2, nu = 3), n = 4)

  expect_equal(noise_moments(nig_noise), noise_nig_moments(mu = 1, sigma = 2, nu = 3))
  expect_equal(noise_moments(gal_noise), noise_gal_moments(mu = 1, sigma = 2, nu = 3))
})

test_that("stationary moment object method rejects non-stationary bases", {
  noise <- noise_nig(
    theta_mu = c(0, 1),
    B_mu = diag(2),
    sigma = 2,
    nu = 3
  )

  expect_error(noise_moments(noise), "stationary")
})
