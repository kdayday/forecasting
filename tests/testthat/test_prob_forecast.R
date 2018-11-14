context("Test forecasting context")

library(forecasting)
library(rvinecopulib)

mock_samp <- function(x) "A sample"
mock_pd <- function(x,y) "A pd"
mock_eval <- function(x,y,z) list(cvar = list(low='low', high='high'), var=list(low=0, high=1))
epsilon=c(0.05, 0.95)

test_that("Basic vine copula forecast initialization is correct.", {
  with_mock(get_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
            vinecop=function(x, ...) NA, marg_transform=function(...) "A transform",
            to_uniform=function(...) "To uniform",
  OUT <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', time=1,  n=3000, epsilon=epsilon))
  expect_identical(OUT$cvar$low, 'low')
  expect_true(is.prob_forecast((OUT)))
  expect_equal(length(OUT), 2)
})

test_that("Marginal distribution optional arguments are passed through.", {
  with_mock(marg_transform=function(x, method='default', anoption=NA) return(list(method=method, anoption=anoption)),
            to_uniform=function(...) return(NA), vinecop=function(...) return(NA), calc_quantiles=mock_pd,
            get_samples= mock_samp, calc_cvar=mock_eval,
            out <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'TX', 'time',
                                  training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Vine copula forecast initialization throws errors", {
  with_mock(get_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
            vinecop=function(x, ...) NA,
            expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', 1,  n='three', epsilon=epsilon)))
  expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=1), 'Odessa', 1,  n=3000, epsilon=epsilon))
})

test_that("CVAR estimate is correct", {
  OUT <- calc_cvar(0:100, c(0.05, 0.95))
  expect_equal(OUT$cvar$low, 2.5)
  expect_equal(OUT$var$low, 5)
  expect_equal(OUT$cvar$high, 97.5)
  expect_equal(OUT$var$high, 95)
})

test_that("CVAR calculation throws errors", {
  expect_error(calc_cvar(matrix(c(0,0,0,0), ncol=2), c(0.5, 2)))
  expect_error(calc_cvar(matrix(c(0,0,0,0), ncol=2), c(0.5, -0.1)))
})

test_that("Interval score calculation is correct.", {
  q <- seq(0, 100, 10)
  names(q) <-sapply(q, FUN=paste, '%', sep='')
  dat <- list(quantiles=q)
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_equal(calc_is(fake_forecast, actual=15, alpha=0.2), 80)
  expect_equal(calc_is(fake_forecast, actual=9, alpha=0.2), 80+2/0.2)
  expect_equal(calc_is(fake_forecast, actual=95, alpha=0.2), 80+2/0.2*5)
})

test_that("Interval score calculation throws error for bad input", {
  q <- seq(0, 100, 10)
  names(q) <-sapply(q, FUN=paste, '%', sep='')
  dat <- list(quantiles=q)
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_error(calc_is(fake_forecast, actual=95, alpha=10), "Alpha.*")
  expect_error(calc_is(fake_forecast, actual=95, alpha=0.1), "Requested quantile is not in the forecast's list of quantiles.") #Unlisted quantile
})

fake_sampling <- function(x,y){
  return(matrix(c(0.5, 0.1, 0.5, 0.2), ncol=2))
}
fake_uniform <- function(x,y){
  return(x*y)
}

test_that('Vine copulas samples are added correctly', {
  dat <- list(n=2, model=NA, d=2, results_transforms=list(5,10))
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  with_mock(rvinecop=fake_sampling, from_uniform=fake_uniform,
            expect_identical(get_samples(fake_forecast), c(7.5, 2.5)))
})
