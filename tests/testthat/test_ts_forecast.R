context("Test forecasting context")

library(forecasting)
library(lubridate)
library(scoringRules)

mock_calc <- function(x, sun_up, start_time, time_step, scale, location, method, n, epsilon) "Calculated"
mock_sunup <- function(x) "Sunup"
fake_class <- function(x,y,z,q, m) sum(x)
fake_class2 <- function(x,y) return(function(x,y,z,q, m) "A forecast")
mock_get <- function(x,y) return(fake_class)
start_time <- ymd(20161130)
time_step <- 0.25
x_site <- list(matrix(c(0,1), ncol = 1))
x_multi <- list(matrix(c(0,1), ncol=2), matrix(c(2,3), ncol=2))
sun_up <- c(TRUE, TRUE)

test_that("ts_forecast object initialization", {
  with_mock(calc_forecasts = mock_calc, check_sunup = mock_sunup,
  OUT <- ts_forecast(x_site, start_time, time_step, 'Site', 'Odessa', 'vine', 3000, 0.05))
  expect_identical(OUT$forecasts, 'Calculated')
  expect_identical(class(OUT), 'ts_forecast')
})

test_that("ts_forecast initialization throws errors", {
  with_mock(calc_forecasts = mock_calc, check_sunup = mock_sunup,
            expect_error(ts_forecast(x_multi, start_time, time_step, 'site', 'Odessa', 'vine', 3000, 0.05)))
  with_mock(calc_forecasts = mock_calc, check_sunup = mock_sunup,
            expect_error(ts_forecast(x_site, start_time, time_step, 'region', 'Odessa', 'vine', 3000, 0.05)))
})

test_that("ts_forecast calculation inserts NA's when sun is down", {
  with_mock(get_forecast_class=fake_class2,
            OUT <- calc_forecasts(x_multi, c(FALSE, TRUE), start_time, time_step, 'region', 'Odessa', 'vine', 3000, 0.05))
  expect_identical(OUT, list(NA, "A forecast"))
})

test_that("ts_forecast class lookup is correct", {
  expect_identical(get_forecast_class('S', 'vine'), prob_1d_site_forecast)
  expect_identical(get_forecast_class('r', 'vine'), prob_nd_vine_forecast)
  expect_identical(get_forecast_class('Total', 'gaussian'), prob_nd_gaussian_forecast)
  expect_identical(get_forecast_class('t', 'E'), prob_nd_empirical_forecast)
  expect_error(get_forecast_class('T', 't'))
})

test_that("ts_forecast calculate list of forecasts", {
  with_mock(get_forecast_class=mock_get,
  OUT <- calc_forecasts(x_multi, sun_up, start_time, time_step, 'site', "Odessa", 'vine', 3000, 0.05))
  expect_identical(OUT, list(1, 5))
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
