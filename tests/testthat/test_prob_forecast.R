context("Test forecasting context")

library(forecasting)
library(VineCopula)

mock_samp <- function(x) "A sample"
mock_pd <- function(x,y) "A pd"
mock_eval <- function(x,y,z) list(cvar = list(low='low', high='high'), var=list(low=0, high=1))
epsilon=c(0.05, 0.95)

test_that("Vine copula forecast is initialized correctly", {
  with_mock(get_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
            pobs=function(x) NA, RVineStructureSelect=function(x, indeptest, familyset) NA,
  OUT <- prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', 1,  3000, epsilon))
  expect_identical(OUT$cvar$low, 'low')
  expect_true(is.prob_forecast((OUT)))
  expect_equal(length(OUT), 2)
})

test_that("Vine copula forecast initialization throws errors", {
  with_mock(get_samples = mock_samp, calc_quantiles=mock_pd, calc_cvar = mock_eval,
            pobs=function(x) NA, RVineStructureSelect=function(x, indeptest, familyset) NA,
            expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=2), 'Odessa', 1,  'three', epsilon)))
  expect_error(prob_nd_vine_forecast(matrix(c(0,0,0,0), ncol=1), 'Odessa', 1,  3000, epsilon))
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

fake_sampling <- function(x,y){
  return(matrix(c(0.5, 0.1, 0.5, 0.1), ncol=2))
}
dat <- list(n=2, model=NA, d=2, training_mat=cbind(0:50, seq(from=0, to=100, by=2)))
fake_forecast <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))

test_that('Vine copulas are sampled and added correctly', {
  with_mock(RVineSim=fake_sampling,
            expect_identical(get_samples(fake_forecast), c(75, 15)))
})
