context("Test forecasting context")

library(forecasting)
library(lubridate)

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

test_that("Integrate telemetry throws error on inputs of wrong length.", {
  expect_error(equalize_telemetry_forecast_length(tel=c(10, 15), fc=c(3, 4, 5)), '*multiples*')
})

test_that("Integrate telemetry calculation is correct", {
  expect_equal(equalize_telemetry_forecast_length(tel=1:4, fc=5:8, agg=TRUE), list(tel=1:4, fc=5:8, tel_2_fc=1))
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, 2, 3, 9, 9), fc=1:2, agg=TRUE), list(tel=c(1, 7), fc=1:2, tel_2_fc=3))
  expect_equal(equalize_telemetry_forecast_length(tel=1:4, fc=5:8, agg=FALSE), list(tel=1:4, fc=5:8, tel_2_fc=1))
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, 2, 3, 9, 9), fc=1:2, agg=FALSE), list(tel=c(1, 0, 2, 3, 9, 9), fc=c(1,1,1,2,2,2), tel_2_fc=3))
})

test_that("Calculation of number of time points without night-time and NaN's in telemetry is correct", {
  fake_x <- structure(list(sun_up=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)), class = c("ts_forecast"))
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(3, 4, 5, NaN, 0, NaN, NaN, 0, 1, -1, NaN, NaN))),
  result <- get_sundown_and_NaN_stats(fake_x, tel=NA))
  expect_equal(result$"Sun-up forecasts", 6)
  expect_equal(result$"Sun-down forecasts", 6)
  expect_equal(result$"Telemetry missing when sun up", 2)
  expect_equal(result$"Telemetry available when sun up", 4)
  expect_equal(result$"Telemetry is 0 when sun forecasted down", 1)
  expect_equal(result$"Telemetry is non-zero when sun forecasted down", 2)
  expect_equal(result$"Telemetry is NaN when sun forecasted down", 3)
})

test_that("Integrate telemetry calculation handles NA's correctly", {
  expect_equal(equalize_telemetry_forecast_length(tel=c(1:4, NaN), fc=1:5), list(tel=c(1:4,NaN), fc=1:5, tel_2_fc=1)) # NA's get passed through
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, NaN, 3, 9, 9, NaN, NaN, NaN), fc=1:3), list(tel=c(NaN, 7, NaN), fc=1:3, tel_2_fc=3)) # Time periods with at least one NA get NAN
})

x1 <- structure(list(quantiles=seq(10, 90, 10), n=2), class="prob_forecast")
names(x1$quantiles) <-sapply(x1$quantiles, FUN=paste, '%', sep='')
x2 <- structure(list(quantiles=seq(5, 45, 5), n=2), class="prob_forecast")
names(x2$quantiles) <-sapply(x1$quantiles, FUN=paste, '%', sep='')
ts <- structure(list(forecasts=list(x1, x1, x2), sun_up=c(FALSE, TRUE, TRUE)), class='ts_forecast')

test_that("Time series of values at certain quantile is extracted correctly.", {
  expect_equal(get_quantile_time_series(ts, 20), c(20, 20, 10))
})

test_that("Extraction of time-series at certain quantile throws error.", {
  expect_error(get_quantile_time_series(ts, 2))
})

test_that("Avg CRPS score handles telemetry longer than forecast", {
  tel <- c(25, 30, 35)
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=tel, tel_2_fc=2)),
            crps=function(x, y) return(y),
                        expect_equal(eval_avg_crps(ts, tel=tel, agg=FALSE), list(sd=stats::sd(c(30, 35)), mean=32.5)),
            expect_equal(eval_avg_crps(ts, tel=tel, agg=TRUE), list(sd=stats::sd(c(30, 35)), mean=32.5)))
})

test_that("Brier score is as expected", {
  with_mock(get_quantile_time_series=function(x,y) return(c(20, 20, 10)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(0, 40, 5), fc=ts)),
    expect_equal(eval_brier(ts, tel=NA, alpha=.8), 0.68))
})

test_that("Brier score calculation handles NaN's", {
  fake_ts <- structure(list(forecasts=list(x1, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(get_quantile_time_series=function(x,y) return(c(20, 20, 10, 10)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(0, 40, 5, NA))),
            expect_equal(eval_brier(fake_ts, tel=NA, alpha=.8), 0.68))
  })

test_that("MAE calculation is correct.", {
  with_mock(get_quantile_time_series=function(x,y) return(c(50, 50, 25)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(60, 60, 5))),
    expect_equal(eval_mae(ts, tel=NA), 15)) # error = 10 and 20
})

test_that("MAE calculation handles NaN's.", {
  fake_ts <- structure(list(forecasts=list(x1, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(get_quantile_time_series=function(x,y) return(c(50, 50, 25, 30)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(60, 60, 5, NaN))),
            expect_equal(eval_mae(fake_ts, tel=NA), 15)) # error = 10 and 20
})

test_that("Average interval score calculation is correct.", {
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(50, 50, 25), tel_2_fc=1)),
            expect_equal(eval_avg_is(ts, tel=NA, alpha=0.2), 60)) # 60=mean(40, 80)
})

test_that("Average interval score calculation handles NaN's.", {
  fake_ts <- structure(list(forecasts=list(x1, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(50, 50, 25, NaN), tel_2_fc=1)),
            expect_equal(eval_avg_is(fake_ts, tel=NA, alpha=0.2), 60)) # 60=mean(80, 40)
})


test_that("Average interval score calculation handles telemetry aggregation (agg=TRUE).", {
  fake_ts <- structure(list(forecasts=list(x1, x2), sun_up=c(TRUE, TRUE)), class='ts_forecast')
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(50, 25), fc=ts, tel_2_fc=2)),
            expect_equal(eval_avg_is(fake_ts, tel=NA, alpha=0.2, agg=TRUE), 60)) # 60=mean(80, 40)
})

test_that("Average interval score calculation handles telemetry longer than forecast (agg=FALSE).", {
  fake_ts <- structure(list(forecasts=list(x1, x2), sun_up=c(TRUE, TRUE)), class='ts_forecast')
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(50, 50, 25, NaN), fc=rep(ts, each=2), tel_2_fc=2)),
            expect_equal(eval_avg_is(fake_ts, tel=NA, alpha=0.2, agg=FALSE), 200/3)) # 200/3=mean(80, 80, 40)
})

test_that("Quantile reliability calculation is correct.", {
  ts <- structure(list(forecasts=list(x1, x1, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE, TRUE)), class='ts_forecast')

  # Hits: quantile 1, 4, 7, 10
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(NaN, 0, 33, 33, 100), fc=ts, tel_2_fc=1)),
    res <- get_quantile_reliability(ts, tel=NA))
  expect_equal(res$quantiles, seq(0.1, 1, by=0.1))
  expect_equal(res$hit_rate, c(0.25, 0, 0, 0.25, 0, 0, 0.25, 0, 0, 0.25))
})

test_that("Quantile reliability calculation handles aggregated telemetry.", {
  ts <- structure(list(forecasts=list(x1, x1, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE, TRUE)), class='ts_forecast')

  # Hits: quantile 1, 4, 7, 10
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(NaN, 0, 33, 33, 100), fc=ts, tel_2_fc=2)),
            res <- get_quantile_reliability(ts, tel=NA, agg=TRUE))
  expect_equal(res$hit_rate, c(0.25, 0, 0, 0.25, 0, 0, 0.25, 0, 0, 0.25))
})

test_that("Quantile reliability calculation handles telemetry longer than forecast (agg=FALSE).", {
  ts <- structure(list(forecasts=list(x1, x1, x2), sun_up=c(FALSE, TRUE, TRUE)), class='ts_forecast')
  # Hits: quantile 2, 4, 3, 7
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(NaN, 0, 12, 33, 12, 33), fc=rep(ts, each=2), tel_2_fc=2)),
            res <- get_quantile_reliability(ts, tel=NA, agg=FALSE))
  expect_equal(res$hit_rate, c(0, 0.25, 0.25, 0.25, 0, 0, 0.25, 0, 0, 0))
})

