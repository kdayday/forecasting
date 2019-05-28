context("Test forecasting context")

library(forecasting)

fake_data <- matrix(c(1, 2, 1, 10, 20, 10), nrow=3, ncol=2)
colnames(fake_data) <- c('irr', 'press')

weights <- c(0.8, 0.2)
h_data <- matrix(c(1, 2.2, 1, 2, 1.1, 11, 20, 10, 22, 10), nrow=5, ncol=2)
colnames(h_data) <- c('irr', 'press')
h_real <- c(1.5, 1, 5, 2, 0)
an1_irr <- weights[1]/sd(h_data[,'irr'])*0.2
an1_press <- weights[2]/sd(h_data[,'press'])*1
an1 <- an1_irr + an1_press
an2 <- weights[1]/sd(h_data[,'irr'])*0.1 + weights[2]/sd(h_data[,'press'])*2
fake_distance <- c(NA, 0.2, 0.3, 0.1, NA)

test_that("Get analogs throws errors. ", {
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 1, c(0.1))) # Unequal weights and features
  expect_error(get_historical_analogs(fake_data[1:2,], NULL, NULL, 1,weights)) # Even number of time points
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 1, c(0.7, 0.2))) # Weights don't add to 1
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 0, weights)) # Fewer than 1 analogs
})

test_that("Analog selection is correct.", {
  with_mock(delle_monache_distance=function(t, ...) return(fake_distance[t]),
            out <- get_historical_analogs(fake_data, h_data, h_real, 2, weights))
  expect_equal(out$obs, c(2, 1))
  expect_equal(out$distance, c(0.1, 0.2))
  expect_equal(out$forecast, h_data[c(4, 2),])
})

test_that("Analog selection handles NA's in metrics.", {
  h_data <- matrix(c(1, NA, 1, 2, 1.1, 11, 20, 10, 22, 10), nrow=5, ncol=2)
  with_mock(delle_monache_distance = function(t, f, h, ...) return(ifelse(all(!is.na(h[t,])), fake_distance[t], NA)),
            out <- get_historical_analogs(fake_data, h_data, h_real, 2, weights))
  expect_equal(out$obs, c(2, 5))
})

test_that("Analog selection only searches time points with available observations.", {
  h_real <- c(1.5, NA, 5, 2, 0)
  with_mock(delle_monache_distance = function(t, f, h, ...) return(ifelse(all(!is.na(h[t,])), fake_distance[t], NA)),
            out <- get_historical_analogs(fake_data, h_data, h_real, 2, weights))
  expect_equal(out$obs, c(2, 5))
})

test_that("Delle Monache feature distance calculation is correct", {
  expect_equal(feature_distance(1, weights[1], sd(h_data[,1]), fake_data, h_data[1:3,]), an1_irr)
})

test_that("Delle Monache feature distance is 0 for matching NAs", {
  h_data <- matrix(c(NA, NA, 1, 2, 1.1, 11, 20, 10, 22, 10), nrow=5, ncol=2)
  fake_data <- matrix(c(NA, NA, 1, 10, 20, 10), nrow=3, ncol=2)
  expect_equal(feature_distance(1, weights[1], sd(h_data[,1], na.rm=T), fake_data, h_data[1:3,]), 0)
})

test_that("Delle Monache feature distance is NA for non-matching NAs", {
  h_data <- matrix(c(NA, NA, 1, 2, 1.1, 11, 20, 10, 22, 10), nrow=5, ncol=2)
  fake_data <- matrix(c(NA, 2.2, 1, 10, 20, 10), nrow=3, ncol=2)
  expect_true(is.na(feature_distance(1, weights[1], sd(h_data[,1], na.rm=T), fake_data, h_data[1:3,])))
})

test_that("Delle Monache total distance calculation is correct", {
  expect_equal(delle_monache_distance(2, fake_data, h_data, weights, c(sd(h_data[,'irr']), sd(h_data[,'press']))), an1)
})

test_that("Delle Monache total distance is missing for edge cases", {
  expect_true(is.na(delle_monache_distance(1, fake_data, h_data, weights, c(sd(h_data[,'irr']), sd(h_data[,'press'])))))
  expect_true(is.na(delle_monache_distance(5, fake_data, h_data, weights, c(sd(h_data[,'irr']), sd(h_data[,'press'])))))
})

test_that("Delle Monache total distance is NA if a feature is NA", {
  h_data <- matrix(c(NA, NA, 1, 2, 1.1, 11, 20, 10, 22, 10), nrow=5, ncol=2)
  fake_data <- matrix(c(NA, 2.2, 1, 10, 20, 10), nrow=3, ncol=2)
  sigmas <- apply(h_data, 2, sd, na.rm=T)
  expect_true(is.na(delle_monache_distance(2, fake_data, h_data, weights, sigmas)))
})

