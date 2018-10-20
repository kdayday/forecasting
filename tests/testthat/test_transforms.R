context("Test forecasting context")

library(forecasting)
library(KernSmooth)


test_that('KDE of CDF is correct', {
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
