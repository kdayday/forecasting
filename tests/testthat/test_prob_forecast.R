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
  expect_true(is.prob_nd_vine_forecast((OUT)))
  expect_equal(length(OUT), 2)
})

test_that("Marginal distribution optional arguments are passed through.", {
  # This actually goes through calc_transforms as well now, which I haven't bothered changing.
  with_mock(marg_transform=function(x, cdf.method='default', anoption=NA) return(list(method=cdf.method, anoption=anoption)),
            to_uniform=function(...) return(NA), vinecop=function(...) return(NA), calc_quantiles=mock_pd,
            get_1d_samples= mock_samp, CVAR=mock_eval,
            out <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'TX', 'time',
                                  training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Marginal distribution optional arguments are passed through and prob forecast optional arguments are extracted.", {
  fake_forecast_class_call <- function(x, location='TX', time='a time', n=3000,  ...){
    prob_nd_vine_forecast(x, location=location, time=time, n=n, ...)}
  with_mock(marg_transform=function(x, cdf.method='default', anoption=NA) return(list(method=cdf.method, anoption=anoption)),
            to_uniform=function(...) return(NA), vinecop=function(...) return(NA), calc_quantiles=mock_pd,
            get_1d_samples= mock_samp, CVAR=mock_eval,
            out <- fake_forecast_class_call(matrix(c(0,0,0,0), ncol=2), location='TX', time='sometime',
                               training_transform_type="geenens", results_transform_type='geenens', anoption=20))
  expect_equal(out$training_transforms[[1]]$anoption, 20)
  expect_equal(out$training_transforms[[1]]$method, 'geenens')
})

test_that("Vine copula forecast initialization throws errors", {
  with_mock(get_1d_samples = mock_samp, calc_quantiles=mock_pd, CVAR = mock_eval,
            vinecop=function(x, ...) NA,
            expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', 1,  n='three')))
  expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=1), 'Odessa', 1,  n=3000))
})

test_that("get_transform_with_unique_xmin_max handles optional arguments correctly", {
  # Test xmin and xmax are lists
  mock_trans <- function(x, cdf.method, ...) {return(list(...))}
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), cdf.method='a method', athing='a', xmin=c(4, 5), xmax=8))
  expect_equal(out, list(athing='a', xmin=5, xmax=8))
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), cdf.method='a method', athing='a', xmin=4, xmax=c(5, 6)))
  expect_equal(out, list(athing='a', xmin=4, xmax=6))
  # No xmin or xmax
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), cdf.method='a method', athing='a'))
  expect_equal(out, list(athing='a'))
  # No arguments
  with_mock(marg_transform=mock_trans,
            out <- get_transform_with_unique_xmin_max(2, matrix(c(0, 1, 2, 3), ncol=2), cdf.method='a method'))
  expect_equal(length(out), 0)
})

test_that("calc_transforms handles inputs correctly", {
  mock_get_trans <- function(idx, dat, cdf.method, ...) {if (cdf.method=='a') return(sum(dat[,idx])) else return(diff(dat[,idx]))}
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
  OUT <- CVAR(fake_x, samples=0:100, epsilon=c(0.05, 0.95))
  expect_equal(OUT, list(cvar=list(low=2.5, high=97.5), var=list(low=5, high=95)))
  # Test default sampling
  with_mock(get_1d_samples=function(...) return(0:100),
            OUT2 <-CVAR(fake_x, epsilon=c(0.05, 0.95)))
  expect_equal(OUT2$var$high, 95)
})

test_that("CVAR calculations throws errors", {
  fake_x <- structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_error(CVAR(fake_x, samples=matrix(c(0,0,0,0), ncol=2), epsilon=c(0.5, 2)), "Bad input*")
  expect_error(CVAR(fake_x, samples=matrix(c(0,0,0,0), ncol=2), epsilon=c(0.5, -0.1)), "Bad input*")

  fake_x <- structure(list(), class = c("prob_forecast", "prob_1d_kde_forecast"))
  expect_error(CVAR(fake_x,epsilon=c(0.5, 2)), "Bad input*")
  expect_error(CVAR(fake_x, epsilon=c(0.5, -0.1)), "Bad input*")
})

test_that("Vine copula quantile calculation adjusts for inputs", {
  fake_x <- structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast"))
  with_mock(quantile=function(x,...) return(sum(x)) , get_1d_samples=function(...) return(1:4),
            OUT <- calc_quantiles(fake_x))
  expect_equal(OUT, 10)
  with_mock(quantile=function(x,...) return(sum(x)) , get_1d_samples=function(...) return(1:4),
            OUT <- calc_quantiles(fake_x, samples=rep(1, times=4), quantiles=seq(0.2, 0.8, by=0.2)))
  expect_equal(OUT, 4)
})

test_that("Vine copula quantile throws error", {
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast")), quantiles=seq(0, 0.75, by=0.25)), "Bad input*")
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_nd_vine_forecast")), quantiles=seq(0.25, 1, by=0.25)), "Bad input*")
})

test_that("Interval score calculation is correct.", {
  q <- seq(0, 100, 10)
  names(q) <-sapply(q, FUN=paste, '%', sep='')
  dat <- list(quantiles=q)
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_equal(IS(fake_forecast, actual=15, alpha=0.2), 80)
  expect_equal(IS(fake_forecast, actual=9, alpha=0.2), 80+2/0.2)
  expect_equal(IS(fake_forecast, actual=95, alpha=0.2), 80+2/0.2*5)
})

test_that("Interval score calculation throws error for bad input", {
  q <- seq(0, 100, 10)
  names(q) <-sapply(q, FUN=paste, '%', sep='')
  dat <- list(quantiles=q)
  fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))
  expect_error(IS(fake_forecast, actual=95, alpha=10), "Alpha.*")
  expect_error(IS(fake_forecast, actual=95, alpha=0.1), "Requested quantile is not in the forecast's list of quantiles.") #Unlisted quantile
})

test_that("CRPS estimation is correct", {
  # Squared: 0.04, 0.16, | 0.36, 0.16
  # 0.04 + 0.5*0.12 = 0.1; 0.16 + 0.5*0.2 + 0.04 * 0.5*0.12= 0.26 + 0.1=0.36
  # 0.1 + 0.36 = 0.46
  with_mock(calc_quantiles=function(...) {return(1:4)},
            out <- CRPS(x=NA, tel=2, quantiles=seq(0.20, 0.8, length.out = 4)))
  expect_equal(out, 0.46)
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

# -----------------------------------------------------------------------------------------------------------
# 1d rank forecast tests

test_that("Quality control function works correctly and throws errors", {
  # Format error
  expect_error(qc_input(matrix(c(2, 4, 3))), "Input data*")
  # Not enough values error
  expect_error(qc_input(c(1, NA, NA)), "*non-NA*")
  expect_equal(qc_input(c(3, NA, 4.5, NA)), c(3, 4.5))
})


test_that("1d rank forecast initialization correctly handles NA's and multiple value.", {
  with_mock(calc_quantiles=mock_pd, qc_input=function(x) return(x),
            OUT <- prob_1d_rank_forecast(c(2, 2, 4, 3, 4), location='Odessa', time=1))
  expect_true(is.prob_1d_rank_forecast((OUT)))
  expect_equal(length(OUT), 1)
  expect_equal(OUT$rank_quantiles$x, c(2, 2, 3, 4, 4))
  expect_equal(OUT$rank_quantiles$y, c(0, 0.25, 0.5, 0.75, 1))
})

test_that('1d rank forecast quantile calculation is correct', {
  fake_forecast <- structure(list(rank_quantiles=list(x=c(1, 6, 11), y=c(0, 0.5, 1))), class = c("prob_forecast", "prob_1d_rank_forecast"))
  OUT <- calc_quantiles(fake_forecast, quantiles=seq(0.25, 0.75, by=0.25))
  expect_equal(unname(OUT), seq(1+2.5, 11-2.5, by=2.5))
  expect_equal(names(OUT), c("25%", "50%", "75%"))
})

test_that('1d rank forecast quantile calculation throws errors', {
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_1d_rank_forecast")), quantiles=seq(0, 0.75, by=0.25)), "Bad input*")
  expect_error(calc_quantiles(structure(list(), class = c("prob_forecast", "prob_1d_rank_forecast")), quantiles=seq(0.25, 1, by=0.25)), "Bad input*")
})

test_that("1D KDE forecast initialization is correct.", {
  with_mock(probempirical=function(dat, anoption= 'b', ...) return(paste(anoption, 'model', sep=' ')),
            calc_quantiles=function(...) return(NA), qc_input=function(x) return(x),
            OUT <- prob_1d_kde_forecast(c(2, 4, 3), location='Odessa', time=1, cdf.method='empirical', anoption='a'))
  expect_true(is.prob_1d_kde_forecast((OUT)))
  expect_equal(length(OUT), 1)
  expect_equal(OUT$model, 'a model')
})

test_that('1d kde forecast quantile calculation is correct', {
  fake_forecast <- structure(list(model=list(x=c(0, 5, 10), u=c(0, 0.5, 1))), class = c("prob_forecast", "prob_1d_kde_forecast"))
  OUT <- calc_quantiles(fake_forecast, quantiles=seq(0.25, 0.75, length.out = 3))
  expect_equal(unname(OUT), seq(2.5, 7.5, by=2.5))
  expect_equal(names(OUT), c("25%", "50%", "75%"))
})

test_that("1d KDE CVAR estimate is correct", {
  fake_x <- structure(list(model=list(x=seq(5, 105, by=10), u=seq(0, 1, by=0.1), d=rep(0.01, times=11))), class = c("prob_forecast", "prob_1d_kde_forecast"))
  OUT <- CVAR(fake_x, epsilon=c(0.2, 0.8))

  # Low side: CVAR = (1/0.2)*trapz from (5,0.05) to (25, 0.25) = (1/0.2)*3 = 15
  expect_equal(OUT, list(cvar=list(low=15, high=95), var=list(low=25, high=85)))
})

test_that("1D BMA forecast discrete-continuous model normalization & summation is correct", {
  PoC <- c(0, 0.1, 0.5)
  mem <- c(1, 5, 11) # Last should be truncated
  w <- c(0.5, 0.1, 0.4)
  xseq <- seq(0.25, 0.75, by=0.25)
  mp <- 10
  fake_x <- structure(list(model=list(A0=PoC, A1=NA, A2=NA, B0=NA, B1=NA, C0=NA, w=w, A_transform=NA, B_transform=NA),
                           max_power=mp, members=mem), class = c("prob_forecast", "prob_1d_bma_forecast"))

  with_mock(get_poc = function(x, A, ...) return(A),
            get_rho = function(x, ...) return(x),
            get_gamma = function(x, ...) return(x),
            dbeta_gamma_rho = function(xseq, g, r) return(xseq*g),
            out <- get_discrete_continuous_model(fake_x, xseq=xseq)) #e.g., (0.25, 0.5, 0.75)*0.5 * (1-0.1)
  expect_equal(out$members$PoC, PoC)
  expect_equal(out$PoC, 0.21)
  expect_equal(out$members$dbeta, matrix(c(0.5*xseq*0.1*(1)/mp, 0.1*xseq*0.5*(1-0.1)/mp, 0.4*xseq*1*(1-0.5)/mp), ncol=3))
  expect_equal(out$xseq, xseq*mp)
  mem_sum <- 0.1*1*0.5 + 0.5*0.9*0.1 + 1*0.5*0.4
  expect_equal(out$dbeta, xseq*mem_sum/mp)
})
