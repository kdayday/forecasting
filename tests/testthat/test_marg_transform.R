context("Test forecasting context")

library(forecasting)
library(pracma)

# ----------------------------------------------------------------------
# Constructor tests

test_that("marg_transform constructor assignments are correct",{
  with_mock(probempirical=function(x, ...) return(list('x'=x-1, 'd'=x+2, 'u'=x+1)),
                        out <- marg_transform(c(1, 2, 3), 'empirical'))
  expect_equal(out$x, c(0,1,2))
  expect_equal(out$d, c(3, 4, 5))
  expect_equal(out$u, c(2, 3, 4))
  expect_equal(out$xmin, 0)
  expect_equal(out$xmax, 2)
  expect_equal(out$method, 'empirical')
})

# ----------------------------------------------------------------------
# Transformation function tests

x <- 0:10
probs <-rep(0.1, length(x))
test_obj <- structure(list('x'=x, 'd'=probs, 'u'=pracma::cumtrapz(x, probs), 'xmin'=min(x), 'xmax'=max(x)), class = c("marg_transform"))

test_that("Input to uniform-to-variable transform must be [0,1]", {
  expect_error(from_uniform(test_obj, c(-0.1, 1.1)))
})

test_that("from_uniform transform is correct", {
  expect_equal(from_uniform(test_obj, c(0.25, NaN, 0.75)), c(2.5, NaN, 7.5))
})

test_that("to_uniform transform is correct", {
  expect_equal(to_uniform(test_obj, c(NaN, 3, 8)), c(NaN, 0.3, 0.8))
})

test_that("KDE uniform transform gives warning if evaluation points are outside bounds of the cdf", {
  expect_warning(result <- to_uniform(test_obj, c(-1, 8))) # Check min side
  expect_warning(result <- to_uniform(test_obj, c(3, 11))) # Check max side
  expect_warning(result <- to_uniform(test_obj, c(-1, 11)))
  expect_equal(result, c(0, 1))
})

test_that("KDE probability density transform gives warning if evaluation points are outside bounds of the pdf.", {
  expect_warning(result <- to_probability(test_obj, c(-1, 11)))
  expect_equal(result, c(0, 0))
})

