context("Test forecasting context")

library(scoringRules)

test_that('emos_model throws error', {
  expect_error(emos_model(tel=4:6, ens=matrix(1:4, ncol=2), max_power=10, par_init=NA), "Must have same*")
})

test_that('emos_model calculation is correct for default coefficients', {
  with_mock(optim = function(par, ...) return(list(par=par)),
            OUT <- emos_model(tel=c(10, 11), ens=matrix(c(2, 3, 5, 7), ncol=2), max_power=10))
  expect_equal(OUT$a, 8)
  expect_equal(OUT$b, c(1, NA))
  expect_equal(OUT$c, 5)
  expect_equal(OUT$d, 1)
})

test_that('emos_model calculation is correct for given coefficients', {
  with_mock(optim = function(par, ...) return(list(par=par)),
            OUT <- emos_model(tel=c(10, 11), ens=matrix(c(2, 3, 5, 7), ncol=2), max_power=10, par_init=list(a=1, b=c(2, 3), c=4, d=5)))
  expect_equal(OUT$a, 1)
  expect_equal(OUT$b, c(2, 3))
  expect_equal(OUT$c, 4)
  expect_equal(OUT$d, 5)
})

test_that('tnorm_crps_obj calculation is correct', {
  ens <- matrix(1:4, ncol=2)
  par <- c(0, 2, 1, 2, 3) # a=1, b=c(2,3) c = 0, d = 2
  max_power <- 10
  obs <- c(17, 18)
  # location = 12, 17; s2 = 3, 8; std_dev = 6, 16
  with_mock(var = function(z) prod(z),
            crps_tnorm = function(obs, location, scale, lower, upper) return(sum(c(location, scale, -obs))),
            expect_equal(tnorm_crps_obj(par, obs, ens, max_power), 16))
})

test_that('tnorm_crps_obj throws error', {
  ens <- matrix(1:4, ncol=2)
  par <- c(0, 2, 1, 2, 3, 1)
  expect_error(tnorm_crps_obj(par, obs=NA, ens, max_power=NA), "Mismatch between*")
})
