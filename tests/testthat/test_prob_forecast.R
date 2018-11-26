context("Test forecasting context")

library(forecasting)
library(rvinecopulib)
library(stats)

mock_samp <- function(x) "A sample"
mock_pd <- function(x,y) "A pd"
mock_eval <- function(x,y,z) list(cvar = list(low='low', high='high'), var=list(low=0, high=1))

test_that("Basic vine copula forecast initialization is correct.", {
  with_mock(get_1d_samples = mock_samp, calc_quantiles=mock_pd,
            vinecop=function(x, ...) NA, calc_transforms=function(...) list('training'='tr', 'results'='res'),
            to_uniform=function(...) "To uniform",
  OUT <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', time=1,  n=3000))
  expect_true(is.prob_forecast((OUT)))
  expect_equal(length(OUT), 2)
})

test_that("Marginal distribution optional arguments are passed through.", {
  # This actually goes through calc_transforms as well now, which I haven't bothered changing.
  with_mock(marg_transform=function(x, method='default', anoption=NA) return(list(method=method, anoption=anoption)),
            to_uniform=function(...) return(NA), vinecop=function(...) return(NA), calc_quantiles=mock_pd,
            get_1d_samples= mock_samp, calc_cvar=mock_eval,
            out <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'TX', 'time',
                                  training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Marginal distribution optional arguments are passed through and prob forecast optional arguments are extracted.", {
  fake_forecast_class_call <- function(x, location='TX', time='a time', n=3000,  ...){
    prob_nd_vine_forecast(x, location=location, time=time, n=n, ...)}
  with_mock(marg_transform=function(x, method='default', anoption=NA) return(list(method=method, anoption=anoption)),
            to_uniform=function(...) return(NA), vinecop=function(...) return(NA), calc_quantiles=mock_pd,
            get_1d_samples= mock_samp, calc_cvar=mock_eval,
            out <- fake_forecast_class_call(matrix(c(0,0,0,0), ncol=2), location='TX', time='sometime',
                               training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Vine copula forecast initialization throws errors", {
  with_mock(get_1d_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
            vinecop=function(x, ...) NA,
            expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', 1,  n='three')))
  expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=1), 'Odessa', 1,  n=3000))
})

test_that("get_transform_with_unique_xmin_max handles optional arguments correctly", {
  # Test xmin and xmax are lists
  mock_trans <- function(x, method=method, ...) {return(list(...))}
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), method='a method', athing='a', xmin=c(4, 5), xmax=8))
  expect_equal(out, list(athing='a', xmin=5, xmax=8))
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), method='a method', athing='a', xmin=4, xmax=c(5, 6)))
  expect_equal(out, list(athing='a', xmin=4, xmax=6))
  # No xmin or xmax
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), method='a method', athing='a'))
  expect_equal(out, list(athing='a'))
  # No arguments
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), method='a method'))
  expect_equal(length(out), 0)
})

test_that("calc_transforms handles inputs correctly", {
  mock_get_trans <- function(idx, dat, method, ...) {if (method=='a') return(sum(dat[,idx])) else return(diff(dat[,idx]))}
  with_mock(get_transform_with_unique_xmin_max=mock_get_trans,
            out <-calc_transforms(matrix(c(5, 2, 7, 1), ncol=2), training_transform_type='a', results_transform_type='b', something='else'))
  expect_equal(out$training, list(7, 8))
  expect_equal(out$results, list(-3, -6))
  with_mock(get_transform_with_unique_xmin_max=mock_get_trans,
            out <-calc_transforms(matrix(c(5, 2, 7, 1), ncol=2), training_transform_type='a', results_transform_type='a'))
  expect_equal(out$results, list(7, 8))
})

test_that("CVAR estimate is correct", {
  # Test sample over-ride
  fake_x <- structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast"))
  OUT <- calc_cvar(fake_x, samples=0:100, epsilon=c(0.05, 0.95))
  expect_equal(OUT$cvar$low, 2.5)
  expect_equal(OUT$var$low, 5)
  expect_equal(OUT$cvar$high, 97.5)
  expect_equal(OUT$var$high, 95)
  # Test default sampling
  with_mock(get_1d_samples=function(...) return(0:100),
            OUT2 <-calc_cvar(fake_x, epsilon=c(0.05, 0.95)))
  expect_equal(OUT2$var$high, 95)
})

test_that("CVAR calculation throws errors", {
  fake_x <- structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_error(calc_cvar(fake_x, samples=matrix(c(0,0,0,0), ncol=2), epsilon=c(0.5, 2)), "Bad input*")
  expect_error(calc_cvar(fake_x, samples=matrix(c(0,0,0,0), ncol=2), epsilon=c(0.5, -0.1)), "Bad input*")
})

test_that("Vine copula quantile calculation adjusts for inputs", {
  fake_x <- structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast"))
  with_mock(quantile=function(x,...) return(sum(x)) , get_1d_samples=function(...) return(1:4),
            OUT <- calc_quantiles(fake_x))
  expect_equal(OUT, 10)
  with_mock(quantile=function(x,...) return(sum(x)) , get_1d_samples=function(...) return(1:4),
            OUT <- calc_quantiles(fake_x, samples=rep(1, times=4), quantile_density=0.2))
  expect_equal(OUT, 4)
})

test_that("Vine copula quantile throws error", {
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast")), quantile_density=1), "Bad input*")
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast")), quantile_density=-0.1), "Bad input*")
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
            expect_identical(get_1d_samples(fake_forecast), c(7.5, 2.5)))
})

test_that('get_variable_domain_grid calc is correct with one dimension length', {
  dat <- list(d=3, results_transforms=list(list('xmin'=1, 'xmax'=3),
                                                          list('xmin'=10, 'xmax'=30),
                                                          list('xmin'=100, 'xmax'=300)) )
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  out<- get_variable_domain_grid(fake_forecast, k=3)
  expect_equal(get_variable_domain_grid(fake_forecast, k=3), expand.grid(c(1, 2, 3), c(10, 20, 30), c(100, 200, 300)))
})

test_that('get_variable_domain_grid calc is correct with unique dimension lengths', {
  dat <- list(d=3, results_transforms=list(list('xmin'=1, 'xmax'=3),
                                           list('xmin'=10, 'xmax'=30),
                                           list('xmin'=100, 'xmax'=300)) )
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_equal(get_variable_domain_grid(fake_forecast, k=c(2, 3, 5)), expand.grid(c(1, 3), c(10, 20, 30), c(100, 150, 200, 250, 300)))
})

test_that('get_variable_domain_grid throws errors on bad input', {
  fake_forecast <- structure(list(d=3), class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_error(get_variable_domain_grid(fake_forecast, k=1), "Bad input*") # Too few samples
  expect_error(get_variable_domain_grid(fake_forecast, k=c(2, 2)), "Bad input*") # incorrect number of inputs
})

test_that('get_joint_density_grid calc is correct', {
  dat <- list(n=3, model=NA, d=2, results_transforms=list(1, 10, 100))
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  with_mock(dvinecop=function(x, c) return(rowSums(x)),
            to_uniform=function(c, x) return(x/10),
            to_probability=function(c,x) return(x/c), # Divide by 'results_transforms' value
            get_variable_domain_grid=function(...) expand.grid(c(1, 2), c(10, 20), c(100, 200)),
            out <- get_joint_density_grid(fake_forecast, 2))
  expect_equal(out$d, c(11.1, 22.4, 24.2, 48.8, 42.2, 84.8, 88.4, 177.6))
  expect_equal(dim(out$grid_points), c(8,3))
})
