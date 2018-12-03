context("Test forecasting context")

library(forecasting)
library(lubridate)
library(scoringRules)

mock_calc <- function(x, sun_up, start_time, time_step, scale, location, method, ...) "Calculated"
fake_class <- function(x, location, time, ...) sum(x)
fake_class2 <- function(x, location, time, ...) "A forecast"
start_time <- ymd(20161130)
time_step <- 0.25
x_site <- array(0:3, dim=c(2, 2))
x_multi <-  array(c(0, 2, 1, 3), dim=c(2, 1, 2))
sun_up <- c(TRUE, TRUE)

test_that("ts_forecast object initialization and sun-up calculation", {
  with_mock(calc_forecasts = mock_calc,
  OUT <- ts_forecast(array(c(0, NA, 3), dim=c(3, 1)), start_time, time_step, 'site', 'Odessa', 'rank'))
  expect_identical(OUT$forecasts, 'Calculated')
  expect_identical(class(OUT), 'ts_forecast')
  expect_equal(OUT$sun_up, c(FALSE, FALSE, TRUE))
})

test_that("ts_forecast initialization throws errors", {
  expect_error(ts_forecast(c(1, 2, 3), start_time, time_step, 'site', 'Odessa', 'rank'), "Bad input*")
  expect_error(ts_forecast(array(1:12, dim=c(2, 2, 3)), start_time, time_step, 'site', 'Odessa', 'vine'), "Data and scale mis-match*")
  expect_error(ts_forecast(array(1:12, dim=c(2, 6)), start_time, time_step, 'region', 'Odessa', 'vine'), "Data and scale mis-match*")
  expect_error(ts_forecast(array(1:12, dim=c(12)), start_time, time_step, 'region', 'Odessa', 'vine'), "Data and scale mis-match*")
  })

test_that("ts_forecast calculation inserts NA's when sun is down", {
  with_mock(get_forecast_class= function(x,y) return(list(class=fake_class2, dim='n')),
            OUT <- calc_forecasts(x_multi, c(FALSE, TRUE), start_time, time_step, 'region', 'Odessa', 'vine'))
  expect_identical(OUT, list(NA, "A forecast"))
})

test_that("ts_forecast class lookup is correct", {
  expect_identical(get_forecast_class('S', 'kde')$class, prob_1d_kde_forecast)
  expect_identical(get_forecast_class('S', 'kde')$dim, '1')
  expect_identical(get_forecast_class('S', 'rank')$class, prob_1d_rank_forecast)
  expect_identical(get_forecast_class('r', 'vine')$class, prob_nd_vine_forecast)
  expect_identical(get_forecast_class('r', 'vine')$dim, 'n')
  expect_identical(get_forecast_class('Total', 'gaussian')$class, prob_nd_gaussian_forecast)
  expect_identical(get_forecast_class('t', 'E')$class, prob_nd_empirical_forecast)
  expect_error(get_forecast_class('T', 't'))
})

test_that("ts_forecast calculate list of forecasts", {
  # n-dimensional
  with_mock(get_forecast_class= function(x,y) return(list(class=fake_class, dim='n')),
  OUT <- calc_forecasts(x_multi, sun_up, start_time, time_step, 'T', "Odessa", 'vine'))
  expect_equal(OUT, list(1, 5))
  # 1-dimensional
  with_mock(get_forecast_class= function(x,y) return(list(class=fake_class, dim='1')),
            OUT <- calc_forecasts(x_site, sun_up, start_time, time_step, 'T', "Odessa", 'vine'))
  expect_equal(OUT, list(2, 4))
})

test_that("ts_forecast calculate passes options through.", {
  opt1 <- "thing"
  opt2 <- 22
  fake_class3 <- function(x, location, time, opt1, opt2) {return(list(thing1=opt1, thing2=opt2))}
  with_mock(get_forecast_class=function(x,y) return(list(class=fake_class3, dim='n')),
            out <- calc_forecasts(array(1:3, dim=c(3,1,1)), sun_up=c(TRUE, TRUE, FALSE), start_time=ymd(20160101), time_step=1, scale='site', location='TX', method='rank', opt1=opt1, opt2=opt2))
  expect_equal(out[[1]]$thing1, opt1)
  expect_equal(out[[1]]$thing2, opt2)
})

# This is bad testing, but mainly to convince myself of the logic.
test_that("Check sun-up works correctly", {
  lst <- list(matrix(c(0,0, 0, 0), ncol=2), matrix(c(2,3), ncol=2), matrix(c(0, 0), ncol=1))
  OUT <- unlist(lapply(lst, check_sunup))
  expect_identical(OUT, c(FALSE, TRUE, FALSE))
})

test_that("Average CRPS calculation throws error on input of wrong length.", {
  act <- c(10, 15)
  fake_ts <- structure(list(forecasts = c(1,2,3)), class = "ts_forecast")
  expect_error(eval_avg_crps(fake_ts, act))
})

x1 <- structure(list(quantiles=seq(0, 100, 10), n=2), class="prob_forecast")
names(x1$quantiles) <-sapply(x1$quantiles, FUN=paste, '%', sep='')
x2 <- structure(list(quantiles=seq(0, 50, 5), n=2), class="prob_forecast")
names(x2$quantiles) <-sapply(x1$quantiles, FUN=paste, '%', sep='')
ts <- structure(list(forecasts=list(x1, x1, x2), sun_up=c(FALSE, TRUE, TRUE)), class='ts_forecast')

test_that("Time series of values at certain quantile is extracted correctly.", {
  expect_equal(get_quantile_time_series(ts, 20), c(20, 20, 10))
})

test_that("Extraction of time-series at certain quantile throws error.", {
  expect_error(get_quantile_time_series(ts, 2))
})

test_that("Brier score is as expected", {
  with_mock(get_quantile_time_series=function(x,y) return(c(20, 20, 10)),
    expect_equal(eval_brier(ts, actuals=c(0, 40, 5), 80), 0.68))
})

test_that("MAE calculation is correct.", {
  with_mock(get_quantile_time_series=function(x,y) return(c(50, 50, 25)),
    expect_equal(eval_mae(ts, actuals=c(60, 60, 5)), 15)) # error = 10 and 20
})

test_that("Average interval score calculation is correct.", {
  expect_equal(eval_avg_is(ts, actuals=c(50, 50, 25), alpha=0.2), 60) # 60=mean(40, 80)
})
