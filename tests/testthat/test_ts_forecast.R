context("Test forecasting context")

library(forecasting)
library(lubridate)

mock_calc <- function(x, sun_up, start_time, time_step, scale, location, method, ...) "Calculated"
fake_class <- function(data.input, location, time, ...) sum(data.input)
fake_class2 <- function(data.input, location, time, ...) "A forecast"
start_time <- ymd(20161130)
time_step <- 0.25
x_site <- array(0:3, dim=c(2, 2))
x_multi <-  array(c(0, 2, 1, 3), dim=c(2, 1, 2))
sun_up <- c(TRUE, TRUE)

identity <- function(i) {i}

test_that("ts_forecast object initialization and sun-up calculation", {
  # Includes testing of the new_ts_forecast subfunction
  with_mock(calc_forecasts = mock_calc,
  OUT <- ts_forecast(array(c(0, NA, 3), dim=c(3, 1)), start_time, time_step, 'site', 'Odessa', 'rank'))
  expect_identical(OUT$forecasts, 'Calculated')
  expect_identical(class(OUT), 'ts_forecast')
  expect_equal(OUT$sun_up, c(FALSE, FALSE, TRUE))
})

test_that("ts_forecast object initialization and sun-up calculation with a pre-calculated list of forecasts", {
  forecasts <- list(NA, structure(list(a=1), class="prob_forecast"), NA)
  OUT <- ts_forecast(forecasts, start_time, time_step, 'site', 'Odessa', 'rank')
  expect_equal(OUT$sun_up, c(FALSE, TRUE, FALSE))
})

test_that("ts_forecast array constructor throw errors", {
  expect_error(ts_forecast(array(1:12, dim=c(2, 2, 3)), start_time, time_step, 'site', 'Odessa', 'vine'), "Data and scale mis-match*")
  expect_error(ts_forecast(array(1:12, dim=c(2, 1, 2, 3)), start_time, time_step, 'region', 'Odessa', 'vine'), "Data and scale mis-match*")
  })

test_that("ts_forecast matrix constructor throws error", {
  expect_error(ts_forecast(matrix(1:12, nrow=2), start_time, time_step, 'region', 'Odessa', 'vine'), "Data and scale mis-match*")
})

test_that("ts_forecast list constructor throws error", {
  expect_error(ts_forecast(list(NA, NA, NA), start_time, time_step, 'region', 'Odessa', 'vine'), "*list of forecasts.")
})

test_that("ts_forecast constructor helper throws error", {
  expect_error(new_ts_forecast(x_multi, start_time, time_step, 'region', 'Odessa', 'vine', sun_up_threshold=0.5, MoreTSArgs = list("time"=2)), "MoreTSArgs*")
})

test_that("ts_forecast calculation inserts NA's when sun is down", {
  with_mock(get_forecast_class= function(x,y) return(list(class=fake_class2, dim='n')),
            OUT <- calc_forecasts(x_multi, c(FALSE, TRUE), start_time, time_step, 'region', 'Odessa', 'vine'))
  expect_identical(OUT, list(NA, "A forecast"))
})

test_that("ts_forecast class lookup is correct", {
  expect_identical(get_forecast_class('Site', 'kde')$class, fc_kde)
  expect_identical(get_forecast_class('Site', 'kde')$dim, '1')
  expect_identical(get_forecast_class('Site', 'rank')$class, fc_binned)
  expect_identical(get_forecast_class('Site', 'bma')$class, fc_bma)
  expect_identical(get_forecast_class('Site', 'climate')$class, fc_climatology)
  expect_error(get_forecast_class('region', 'vine'), "Not implemented*")
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
  fake_class3 <- function(data.input, location, time, opt1, opt2) {return(list(thing1=opt1, thing2=opt2))}
  with_mock(get_forecast_class=function(x,y) return(list(class=fake_class3, dim='n')),
            out <- calc_forecasts(array(1:3, dim=c(3,1,1)), sun_up=c(TRUE, TRUE, FALSE), start_time=ymd(20160101), time_step=1, scale='site', location='TX', method='rank', opt1=opt1, opt2=opt2))
  expect_equal(out[[1]]$thing1, opt1)
  expect_equal(out[[1]]$thing2, opt2)
})

test_that("ts_forecast calculate handles additional time-series arguments correctly ", {
  fake_class3 <- function(data.input, location, time, opt1, opt2) {return(list(thing1=opt1, thing2=opt2))}
  with_mock(get_forecast_class=function(x,y) return(list(class=fake_class3, dim='n')),
            out <- calc_forecasts(array(1:3, dim=c(3,1,1)), sun_up=c(TRUE, TRUE, FALSE), start_time=ymd(20160101), time_step=1, scale='site', location='TX', method='rank',
                                  MoreTSArgs = list(opt1=3:5, opt2=7:9)))
  expect_equal(out[[1]]$thing1, 3)
  expect_equal(out[[2]]$thing2, 8)
})

# This is bad testing, but mainly to convince myself of the logic.
test_that("Check sun-up works correctly", {
  lst <- list(matrix(c(0,0, 0, 0), ncol=2), matrix(c(2,3), ncol=2), matrix(c(0, 0.5), ncol=1), matrix(c(0, 0.4), ncol=1))
  OUT <- unlist(lapply(lst, check_sunup, sun_up_threshold=0.5))
  expect_identical(OUT, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("Equalize telemetry throws error on input checks.", {
  expect_error(equalize_telemetry_forecast_length(tel=c(10, 15), fc=c(3, 4, 5)), '*even multiple*')
  expect_error(equalize_telemetry_forecast_length(tel=1:15, fc=c(3, 4, 5)), '*even multiple*')
  expect_error(equalize_telemetry_forecast_length(tel=4, fc=1:3), '*even multiple*') # Test single value rather than vector inputs
  expect_error(equalize_telemetry_forecast_length(tel=1:4, fc=3), '*even multiple*') # Test single value rather than vector inputs
  expect_error(equalize_telemetry_forecast_length(tel=c(10, 15), fc=c(3, 4), align="end of hour"), 'Unknown method*')
})

test_that("Equalize telemetry does nothing for equal length forecasts", {
  expect_equal(equalize_telemetry_forecast_length(tel=1:4, fc=5:8)[c('tel', 'fc', 'tel_2_fc')],
               list(tel=1:4, fc=5:8, tel_2_fc=1))
  expect_equal(equalize_telemetry_forecast_length(tel=1, fc=5)[c('tel', 'fc', 'tel_2_fc')],
               list(tel=1, fc=5, tel_2_fc=1))
})

test_that("Equalize telemetry calculation is correct for top of hour", {
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, 2, 1, 3, 3, 9, 9), fc=1:2, agg=TRUE, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(1, 6), fc=1:2, tel_2_fc=4))
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, 2, 1, 3, 3, 9, 9), fc=1:2, agg=FALSE, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(1, 0, 2, 1, 3, 3, 9, 9), fc=c(1,1,1,1,2,2,2,2), tel_2_fc=4))
})

test_that("Equalize telemetry calculation handles length rather than vector inputs", {
  expect_equal(equalize_telemetry_forecast_length(tel=c(1, 0, 2, 1, 3, 3, 9, 9), fc=2, agg=TRUE, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(1, 6), fc=2, tel_2_fc=4))
  expect_equal(equalize_telemetry_forecast_length(tel=8, fc=1:2, agg=FALSE, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=8, fc=c(1,1,1,1,2,2,2,2), tel_2_fc=4))
})

test_that("Equalize telemetry calculation is correct for half-hour alignment", {
  expect_equal(equalize_telemetry_forecast_length(tel=1:12, fc=1:3, agg=TRUE, align="half-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(NaN, 18, 34)/4, fc=1:3, tel_2_fc=4))
  expect_equal(equalize_telemetry_forecast_length(tel=1:12, fc=1:3, agg=FALSE, align="half-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=1:12, fc=c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, NA, NA), tel_2_fc=4))
})

test_that("Equalize telemetry output translation functions are correct", {
  func <- equalize_telemetry_forecast_length(tel=1:8, fc=1:2, agg=TRUE, align="end-of-hour")$translate_forecast_index
  expect_equal(func(1:12), 1:12)

  func <- equalize_telemetry_forecast_length(tel=1:8, fc=1:2, agg=FALSE, align="end-of-hour")$translate_forecast_index
  expect_equal(func(1:12), c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))

  func <- equalize_telemetry_forecast_length(tel=1:8, fc=1:2, agg=TRUE, align="half-hour")$translate_forecast_index
  expect_equal(func(1:12), 1:12)

  func <- equalize_telemetry_forecast_length(tel=1:8, fc=1:2, agg=FALSE, align="half-hour")$translate_forecast_index
  expect_equal(func(1:12), c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4))
})

test_that("Equalize telemetry output translation functions are correct for length rather than vector inputs", {
  func <- equalize_telemetry_forecast_length(tel=1:8, fc=2, agg=TRUE, align="end-of-hour")$translate_forecast_index
  expect_equal(func(1:12), 1:12)

  func <- equalize_telemetry_forecast_length(tel=8, fc=1:2, agg=FALSE, align="end-of-hour")$translate_forecast_index
  expect_equal(func(1:12), c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))

  func <- equalize_telemetry_forecast_length(tel=1:8, fc=2, agg=TRUE, align="half-hour")$translate_forecast_index
  expect_equal(func(1:12), 1:12)

  func <- equalize_telemetry_forecast_length(tel=8, fc=1:2, agg=FALSE, align="half-hour")$translate_forecast_index
  expect_equal(func(1:12), c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4))
})


test_that("Equalize telemetry calculation handles NA's correctly", {
  expect_equal(equalize_telemetry_forecast_length(tel=c(1:4, NaN), fc=1:5, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(1:4,NaN), fc=1:5, tel_2_fc=1)) # NA's get passed through
  expect_equal(equalize_telemetry_forecast_length(tel=c(2, 1, 0, NaN, 3, 3, 9, 9, NaN, NaN, NaN, NaN), fc=1:3, align="end-of-hour")[c('tel', 'fc', 'tel_2_fc')],
               list(tel=c(NaN, 6, NaN), fc=1:3, tel_2_fc=4)) # Time periods with at least one NA get NAN
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

quants <- seq(0.1, 0.9, by=0.1)
x1 <- structure(list(quantiles=list(x=seq(10, 90, 10), q=quants), n=2), class="prob_forecast")
x2 <- structure(list(quantiles=list(x=seq(5, 45, 5), q=quants), n=2), class="prob_forecast")
ts <- structure(list(forecasts=list(NA, x1, x2), sun_up=c(FALSE, TRUE, TRUE)), class='ts_forecast')

test_that("Plot ts_forecast throws error", {
  expect_error(plot(structure(list(NA), class="prob_forecast"), window=c(1,4,5)))
})

test_that("Plot and write index subfunction is correct", {
  expect_equal(get_index_window(21:30, window=NA), 1:10)
  expect_equal(get_index_window(21:30, window=c(2, 8)), 2:8)
})


test_that("Time series of values at certain quantile is extracted correctly.", {
  expect_equal(get_quantile_time_series(ts, 20), c(NA, 20, 10))
})

test_that("Extraction of time-series at certain quantile throws error.", {
  expect_error(get_quantile_time_series(ts, 2))
})

test_that("Equalize normalization factor is correct", {
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=tel, fc=ts)),
    expect_equal(equalize_normalization_factor_length(3, 1:4), c(3,3,3,3)))
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=tel, fc=ts)),
            expect_equal(equalize_normalization_factor_length(5:8, 1:4), 5:8))
})

test_that("Avg CRPS score handles telemetry longer than forecast", {
  tel <- c(25, 30, 35)
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=tel, tel_2_fc=2, translate_forecast_index=identity)),
            CRPS=function(x, tel, ...) return(tel),
            expect_equal(CRPS_avg(ts, tel=tel, agg=FALSE), list(mean=32.5, min=30, max=35, sd=stats::sd(c(30, 35)), median=32.5)),
            expect_equal(CRPS_avg(ts, tel=tel, agg=TRUE), list(mean=32.5, min=30, max=35, sd=stats::sd(c(30, 35)), median=32.5)))
})

test_that("Reliability index is correct", {
  with_mock(calc_PIT_histogram=function(ts, tel, nbins, ...) return(list(bin_hits=c(0, 3, 9))),
            expect_equal(RI(ts=NA, tel=NA, nbins=3), 10/12))
})

test_that("Percentile reliability index is correct", {
  with_mock(calc_PIT_histogram=function(ts, tel, ...) return(list(bin_hits=c(4, 5, 6, 7, 3, rep(5, times=95)))),
            expect_equal(PRI(ts=NA, tel=NA), 0.008))
})

test_that("Brier score is as expected", {
  ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(get_quantile_time_series=function(x,y) return(c(20, 20, 10, 10)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(0, 40, 5, 15), fc=ts, translate_forecast_index=identity)),
    expect_equal(Brier(ts, tel=NA, PoE=.8), 0.24)) # mean( (-0.2)^2, (0.8)^2, (-0.2)^2) = (0.64 +0.04 + 0.04)/3 = 0.24
})

test_that("Brier score calculation handles NaN's", {
  fake_ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(get_quantile_time_series=function(x,y) return(c(20, 20, 10, 10)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(0, 40, 5, NA), translate_forecast_index=identity)),
            expect_equal(Brier(fake_ts, tel=NA, PoE=.8), 0.34))
  })

test_that("Brier score for power threshold is correct", {
  ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  expect_equal(Brier_power(ts, tel=c(0, 40, 5, NA), thresholds=10), 0.325) # mean( (0.1)^2 , (-0.8)^2) = (0.01 + 0.64)/2 = 0.325
})

test_that("Brier score for power handles a vector of thresholds", {
  ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  expect_equal(Brier_power(ts, tel=c(0, 40, 5, NA), thresholds=c(10, 45, 20)), c(0.325, 0.15625, 0.2)) # mean( (-0.55)^2 , (-0.1)^2) = (0.3025 + 0.01)/2 = 0.15625
})

test_that("MAE calculation is correct.", {
  with_mock(get_quantile_time_series=function(x,y) return(c(50, 50, 25)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(60, 60, 5), translate_forecast_index=identity)),
    expect_equal(MAE(ts, tel=NA), 15)) # error = 10 and 20
})

test_that("MAE calculation handles NaN's.", {
  fake_ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(FALSE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(get_quantile_time_series=function(x,y) return(c(50, 50, 25, 30)),
            equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(60, 60, 5, NaN), translate_forecast_index=identity)),
            expect_equal(MAE(fake_ts, tel=NA), 15)) # error = 10 and 20
})

test_that("Average interval score calculation is correct.", {
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(50, 50, 25), tel_2_fc=1, translate_forecast_index=identity)),
            expect_equal(IS_avg(ts, tel=NA, alpha=0.2), list(mean=60, min=40, max=80, sd=stats::sd(c(40, 80))))) # 60=mean(40, 80)
})

test_that("Average interval score calculation handles NaN's.", {
  fake_ts <- structure(list(forecasts=list(NA, x1, x2, x2), sun_up=c(TRUE, TRUE, TRUE, TRUE)), class='ts_forecast')
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(fc=ts, tel=c(50, 50, 25, NaN), tel_2_fc=1, translate_forecast_index=identity)),
            expect_equal(IS_avg(fake_ts, tel=NA, alpha=0.2), list(mean=60, min=40, max=80, sd=stats::sd(c(40, 80))))) # 60=mean(80, 40)
})

test_that("Average interval score calculation handles telemetry longer than forecast (index translation function).", {
  fake_ts <- structure(list(forecasts=list(x1, x2), sun_up=c(TRUE, TRUE)), class='ts_forecast')
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(50, 50, 25, NaN), fc=rep(ts, each=2), tel_2_fc=2,
                                                                                  translate_forecast_index=function(i) {floor((i-1)/2)+1})),
            expect_equal(IS_avg(fake_ts, tel=NA, alpha=0.2, agg=FALSE), list(mean=200/3, min=40, max=80, sd=stats::sd(c(40, 80, 80))))) # 200/3=mean(80, 80, 40)
})

test_that("Average sharpness score calculation is correct.", {
  expect_equal(sharpness_avg(ts, tel=1:3, alpha=0.2, normalize.by = c(2, 10, 10)), list(mean=6, min=4, max=8, sd=stats::sd(c(8,4)))) # 6=mean(8, 4)
})

test_that("Average sharpness score skips times when telemetry is missing.", {
  expect_equal(sharpness_avg(ts, tel=c(1, NA, 2), alpha=0.2, normalize.by = c(2, 20, 10)), list(mean=4, min=4, max=4, sd=stats::sd(c(4)))) # 6=mean(80, 40)
})

test_that("Quantile reliability calculation is correct.", {
  ts <- structure(list(forecasts=list(NA, x1, x1, x2, x2), sun_up=c(TRUE, TRUE, TRUE, TRUE, TRUE)), class='ts_forecast')

  # Hits: quantile 1, 4, 7, 10
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(NaN, 0, 33, 33, 100), fc=ts, tel_2_fc=1, translate_forecast_index=identity)),
    res <- get_quantile_reliability(ts, tel=NA))
  expect_equal(res$quantiles, seq(0.1, 1, by=0.1))
  expect_equal(res$hit_rate, c(0.25, 0, 0, 0.25, 0, 0, 0.25, 0, 0, 0.25))
})

test_that("Quantile reliability calculation handles telemetry longer than forecast (index translation function).", {
  ts <- structure(list(forecasts=list(NA, x1, x2), sun_up=c(TRUE, TRUE, TRUE)), class='ts_forecast')
  # Hits: quantile 2, 4, 3, 7
  with_mock(equalize_telemetry_forecast_length=function(tel, ts, ...) return(list(tel=c(NaN, 0, 12, 33, 12, 33), fc=rep(ts, each=2), tel_2_fc=2,
                                                                                  translate_forecast_index=function(i) {floor((i-1)/2)+1})),
            res <- get_quantile_reliability(ts, tel=NA, agg=FALSE))
  expect_equal(res$hit_rate, c(0, 0.25, 0.25, 0.25, 0, 0, 0.25, 0, 0, 0))
})

