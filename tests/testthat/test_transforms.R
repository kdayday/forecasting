context("Test forecasting context")

library(forecasting)
library(KernSmooth)


test_that('KDE of CDF data assignments are correct', {
  probs <-rep(0.9, 10)
  with_mock(bkde=function(...) return(list('x'=1:10, 'y'=probs)),
            out <- cdf_kde(c(1), bandwidth=1))
  expect_equal(out$x, 1:10)
  expect_equal(out$pdf, probs)
  expect_equal(out$cdf, seq(0.1, 1, by=0.1))
})


test_that("Get transform throws error.", {
  expect_error(get_transform(1:10, transform_method='norm'))
})

test_that("Get transform lookup is correct.", {
  with_mock(cdf_kde=function(...) return('cdf_kde'),
    expect_equal(get_transform(1:10), 'cdf_kde'))
  expect_warning(  with_mock(cdf_em=function(x) return(x),
                             OUT <- get_transform(1:10, transform_method='empirical')))
  expect_equal(OUT, 1:10)
})

probs <-rep(0.9, 10)
with_mock(bkde=function(...) return(list('x'=1:10, 'y'=probs)),
          out <- cdf_kde(c(1), bandwidth=1))

test_that("KDE uniform transform is correct", {
  expect_equal(to_uniform(out, c(3, 8)), c(0.3, 0.8))
})

test_that("KDE uniform transform gives warning if evaluation points are outside bounds of the cdf", {
  expect_warning(result <- to_uniform(out, c(-1, 8))) # Check min side
  expect_warning(result <- to_uniform(out, c(3, 11))) # Check max side
  expect_warning(result <- to_uniform(out, c(-1, 11)))
  expect_equal(result, c(0, 1))
})

test_that("KDE probability density transform gives warning if evaluation points are outside bounds of the pdf.", {
  expect_warning(result <- to_probability(out, c(-1, 11)))
  expect_equal(result, c(0, 0))
})


test_that("Input to uniform-to-variable transform must be [0,1]", {
  expect_error(from_uniform(result, c(-0.1, 1.1)))
})
