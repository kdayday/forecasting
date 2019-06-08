context("Test forecasting context")

library(forecasting)

fake_data <- matrix(c(1.1, 2, 1, 10, 20, 10), nrow=3, ncol=2)

weights <- c(0.8, 0.2)
h_data <- array(c(1, 2.2, 1, 2.2, 1, 2, 1, 2, 1.1,
                  11, 20, 10, 20, 10, 22, 10, 22, 10), dim=c(3, 3, 2))
h_real <- c(1, 5, 2)
an1_irr <- weights[1]/sd(h_data[,,1])*sqrt(0.2^2 + 0.1^2)
an1_press <- weights[2]/sd(h_data[,,2])*1
an1 <- an1_irr + an1_press
an2 <- weights[1]/sd(h_data[,,1])*0.1 + weights[2]/sd(h_data[,,2])*2
fake_distance <- c(0.2, 0.3, 0.1)

test_that("Get analogs throws errors. ", {
  expect_error(get_historical_analogs(fake_data, h_data, h_real, 1, c(0.1)), "Must have*") # Unequal weights and features
  expect_error(get_historical_analogs(fake_data[1:2,], h_data, h_real, 1, weights), "*odd number*") # Even number of time points
  expect_error(get_historical_analogs(fake_data, h_data, h_real, 1, c(0.7, 0.2)), "Weights must sum*") # Weights don't add to 1
  expect_error(get_historical_analogs(fake_data, h_data, h_real, 0, weights), "At least 1*") # Fewer than 1 analogs
  expect_error(get_historical_analogs(fake_data, h_data, c(1.5, 1, 5, 2, 0, 2), 1, weights), "*time resolution.") # Forecast and telemetry of different lengths
})

test_that("Analog selection is correct.", {
  with_mock(delle_monache_distance=function(h, f, ...) return(abs(sum(fake_data - h))),
            out <- get_historical_analogs(fake_data, h_data, h_real, 2, weights))
  expect_equal(out$obs, c(1, 2))
  expect_equal(out$distance, c(1.1, 2))
  expect_equal(out$forecast, h_data[c(1,3),,])
})

test_that("Analog selection handles NA's in metrics.", {
  h_data[1,1,1] <- NA
  with_mock(delle_monache_distance = function(h, f, ...) return(ifelse(all(!is.na(h)), abs(sum(fake_data - h)), NA)),
            out <- get_historical_analogs(fake_data, h_data, h_real, 2, weights))
  expect_equal(out$obs, c(2, 5))
})

test_that("Analog selection only searches time points with available observations & fills insufficient analogs with NaNs", {
  h_real[1:2] <- NA
  with_mock(delle_monache_distance = function(h, f, weights, ...) return(ifelse(all(!is.na(h)), abs(sum(fake_data - h)), NA)),
            out <- get_historical_analogs(fake_data, h_data, h_real, 4, weights))
  expect_equal(out$obs, c(2, NA, NaN, NaN))
  expect_equal(out$distance, c(2, NA, NaN, NaN))
  expect_equal(out$forecast, aperm(array(c(h_data[3,,], rep(NaN, times=18)), dim=c(3, 2, 4)), c(3, 1,2)))
})

test_that("Get analogs throws errors if no valid analogs are found. ", {
  h_real <- c(NA, NA, NA)
  with_mock(delle_monache_distance=function(h, f, ...) return(abs(sum(fake_data - h))),
            expect_error(get_historical_analogs(fake_data, h_data, h_real, 2, weights), "No viable*"))
})

test_that("Delle Monache feature distance calculation is correct", {
  expect_equal(feature_distance(weights[1], sd(h_data[,,1]), fake_data[,1], h_data[1,,1]), an1_irr)
})

test_that("Delle Monache feature distance calculation is correct for time window of length 1", {
  expect_equal(feature_distance(weights[1], sd(h_data[,,1]), fake_data[2,1], h_data[2,1,1]), weights[1]/sd(h_data[,,1])*0.2)
})

test_that("Delle Monache feature distance is 0 for matching NAs", {
  h_data[1, 1:2, 1] <- NA
  fake_data[1:2,1] <- NA
  expect_equal(feature_distance(weights[1], sd(h_data[,,1], na.rm=T), fake_data[,1], h_data[1,,1]), 0)
})

test_that("Delle Monache feature distance is NA for non-matching NAs", {
  h_data[1, 1:2, 1] <- NA
  fake_data[1,1] <- NA
  expect_true(is.na(feature_distance(weights[1], sd(h_data[,,1], na.rm=T), fake_data[,1], h_data[,,1])))
})

test_that("Delle Monache total distance calculation is correct", {
  expect_equal(delle_monache_distance(fake_data, h_data[1,,], weights, c(sd(h_data[,,1]), sd(h_data[,,2]))), an1)
})

test_that("Delle Monache total distance is NA if a feature is NA", {
  h_data[1, 1:2, 1] <- NA
  fake_data[1,1] <- NA
  sigmas <- apply(h_data, 3, sd, na.rm=T)
  expect_true(is.na(delle_monache_distance(fake_data, h_data[1,,], weights, sigmas)))
})

