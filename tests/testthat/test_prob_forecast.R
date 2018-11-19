context("Test forecasting context")

library(forecasting)
library(rvinecopulib)

mock_samp <- function(x) "A sample"
mock_pd <- function(x,y) "A pd"
mock_eval <- function(x,y,z) list(cvar = list(low='low', high='high'), var=list(low=0, high=1))
epsilon=c(0.05, 0.95)

test_that("Basic vine copula forecast initialization is correct.", {
  with_mock(get_1d_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
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
            get_1d_samples= mock_samp, calc_cvar=mock_eval,
            out <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'TX', 'time',
                                  training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Marginal distribution optional arguments are passed through and prob forecast optional arguments are extracted.", {
  fake_forecast_class_call <- function(x, location='TX', time='a time', n=3000, epsilon=c(0.05,0.95), ...){
    prob_nd_vine_forecast(x, location=location, time=time, n=n, epsilon=epsilon, ...)}
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
