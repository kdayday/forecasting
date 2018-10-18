context("Test forecasting context")

library(forecasting)

fake_data <- matrix(c(1, 2, 1, 10, 20, 10, 100, 200, 100), nrow=3, ncol=3)
colnames(fake_data) <- c('irr', 'press', 'precip')

weights <- c(0.8, 0.2)
h_data <- matrix(c(1, 2.2, 1, 2, 1.1, 11, 20, 10, 22, 10, 110, 220, 110, 180, 110), nrow=5, ncol=3)
colnames(h_data) <- c('irr', 'press', 'precip')
h_real <- c(1.5, 1, 5, 2, 0)
an1_irr <- weights[1]/sd(h_data[,'irr'])*0.2
an1_press <- weights[2]/sd(h_data[,'press'])*1
an1 <- an1_irr + an1_press
an2 <- weights[1]/sd(h_data[,'irr'])*0.1 + weights[2]/sd(h_data[,'press'])*2
features <- c('irr', 'press')

test_that("Get analogs throws errors. ", {
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 1, c('feat1', 'feat2'), c(0.1))) # Unequal weights and features
  expect_error(get_historical_analogs(fake_data[1:2,], NULL, NULL, 1, c('feat1', 'feat2'), weights)) # Even number of time points
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 1, c('feat1', 'feat2'), c(0.7, 0.2))) # Weights don't add to 1
  expect_error(get_historical_analogs(fake_data, NULL, NULL, 0, c('feat1', 'feat2'), weights)) # Fewer than 1 analogs
})

test_that("Analog selection is correct.", {
  out <- get_historical_analogs(fake_data, h_data, h_real, 2, features, weights)
  expect_equal(out$real, c(2, 1))
  expect_equal(out$distance, c(an2, an1))
  expect_equal(out$forecast, h_data[c(4, 2),])
})

test_that("Delle Monache feature distance calculation is correct", {
  expect_equal(feature_distance('irr', weights[1], sd(h_data[,'irr']), fake_data, h_data[1:3,]), an1_irr)
})

test_that("Delle Monache total distance calculation is correct", {
  expect_equal(delle_monache_distance(2, fake_data, h_data,  features, weights, c(sd(h_data[,'irr']), sd(h_data[,'press']))), an1)
  expect_equal(delle_monache_distance(4, fake_data, h_data,  features, weights, c(sd(h_data[,'irr']), sd(h_data[,'press']))), an2)
  expect_true(is.na(delle_monache_distance(1, fake_data, h_data,  features, weights, c(sd(h_data[,'irr']), sd(h_data[,'press'])))))
  expect_true(is.na(delle_monache_distance(5, fake_data, h_data,  features, weights, c(sd(h_data[,'irr']), sd(h_data[,'press'])))))
})
