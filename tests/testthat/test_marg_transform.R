context("Test forecasting context")

library(forecasting)
library(KernSmooth)
library(kde1d)

# ----------------------------------------------------------------------
# Constructor tests
test_that("marg_transform constructor throws errors", {
  expect_error(marg_tranform(c(1, 2, 3, 4), method='probtranskde'))
})

test_that("marg_transform constructor assignments are correct",{
  with_mock(probempirical=function(x, ...) return(list('x'=x+2, 'd'=x-1, 'u'=x+1)),
                        out <- marg_transform(c(1, 2, 3), 'empirical'))
  expect_equal(out$x, c(3, 4, 5))
  expect_equal(out$d, c(0,1,2))
  expect_equal(out$u, c(2, 3, 4))
  expect_equal(out$xmin, 3)
  expect_equal(out$xmax, 5)
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
  expect_equal(from_uniform(test_obj, c(0.25, 0.75)), c(2.5, 7.5))
})

test_that("to_uniform transform is correct", {
  expect_equal(to_uniform(test_obj, c(3, 8)), c(0.3, 0.8))
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

# -----------------------------------------------------------------------------
# Test containers for existing distribution estimation functions

test_that("probempirical handles both a given xmax", {
  out <- probempirical(c(1,2, 4), xmax=5)
  expect_equal(out$x, c(0, 1, 2, 4, 5))
  expect_equal(out$u, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("probempirical replaces given xmax with actual data maximum", {
  out <- probempirical(c(1, 6), xmax=5)
  expect_equal(out$x, c(0, 1,  6))
  expect_equal(out$u, c(0, 0.5, 1))
})

test_that("probempirical handles default max", {
  out <- probempirical(c(1, 2, 4, 5))
  expect_equal(out$x, c(0, 1,  2, 4, 5))
  expect_equal(out$u, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("probempirical handles both duplicate values. ", {
  out <- probempirical(c(1, 2, 2, 4), xmax=5)
  expect_equal(out$x, c(0, 1, 2, 4, 5))
  expect_equal(out$u, c(0, 0.2, 0.6, 0.8, 1))
})

test_that("probkde1d container is correct", {
  with_mock(kde1d=function(x, y=NaN) return(list('grid_points'=x, 'values'=y)),
              out <- probkde1d(0:10, y=rep(0.1, 11)))
  expect_equal(out$x, 0:10)
  expect_equal(out$u, seq(0, 1, by=0.1))
})

test_that("probkde container is correct", {
  with_mock(bkde=function(x, y=NaN) return(list('x'=x, 'y'=y)),
            out <- probkde(0:10, y=rep(0.1, 11)))
  expect_equal(out$x, 0:10)
  expect_equal(out$u, seq(0, 1, by=0.1))
})

# -----------------------------------------------------------------------------
# Test auxillary functions for Geenens et. al transform

test_that("probtranskde throws errors", {
  expect_error(probtranskde(x, xmin=NaN, xmax=1), "Bad input.*")
  expect_error(probtranskde(x, xmin=2, xmax=1), "Bad input.*")
})

test_that("transform_to_real_line throws errors", {
  # Negative value errors
  expect_error(transform_to_real_line(c(-1, 2, 4), xmin=0, xmax=NaN, scale=0.5, zero_offset = 0.1), "Transformation*")
  expect_error(transform_to_real_line(c(-1, 2, 4), xmin=0, xmax=10, scale=0.5, zero_offset = 0.1), "Transformation*")
  # zero errors
  expect_error(transform_to_real_line(c(0, 2, 4), xmin=0, xmax=4, scale=0.9, zero_offset = 0), "Transformation*")
  expect_error(transform_to_real_line(c(0, 2, 4), xmin=0, xmax=NaN, scale=0.9, zero_offset = 0), "Transformation*")
  # 1 errors
  expect_error(transform_to_real_line(c(0, 2, 4), xmin=0, xmax=4, scale=1, zero_offset = 0.1), "Transformation*")
})

test_that("get_output_seq is correct.", {
  expect_equal(get_output_seq(x=c(1, 2, 3, 5), xmin=1, xmax=NaN, n.res=4, scaler=2), c(0.2, 0.4, 0.6, 0.8)*8)
  expect_equal(get_output_seq(x=c(1, 2, 3, 5), xmin=0, xmax=3, n.res=4, scaler=2), c(0.2, 0.4, 0.6, 0.8))
})


test_that("Scale_01 adjusts zeros and xmin for log transform.", {
  expect_equal(scale_01(x=c(1, 2, 4), xmin=1, xmax=NaN, scale=0.5, zero_offset = 0.0001), c(0.0001, 1, 3))
  expect_equal(scale_01(x=c(1, 2, 4), xmin=0, xmax=NaN, scale=0.5, zero_offset = 0.0001), c(1, 2, 4))
})

test_that("Scale_01 rescales for probit transform.", {
  expect_equal(scale_01(x=c(1, 2, 4), xmin=1, xmax=10, scale=0.5, zero_offset = 0.0001), c(0.0001, 0.05, 0.15))
})

test_that("Scale functions throw errors for scaling factors not in (0,1]", {
  expect_error(scale_full(x01=c(0, 1, 2), d01=c(2, 3, 4), xmin=0, xmax=NaN, scale=10), "Scale*")
  expect_error(scale_full(x01=c(0, 1, 2), d01=c(2, 3, 4), xmin=0, xmax=NaN, scale=0), "Scale*")
} )

test_that("Scale_full adjusts log transform for xmin.", {
  out <- scale_full(x01=c(0, 1, 2), d01=c(2, 3, 4), xmin=1, xmax=NaN, scale=0.9)
  expect_equal(out$x, c(1, 2, 3))
  expect_equal(out$d, c(2, 3, 4))
})

test_that("Scale_full rescales for probit transform.", {
  out <- scale_full(x01=c(0, 0.25, 0.5), d01=c(2, 3, 4), xmin=1, xmax=10, scale=0.5)
  expect_equal(out$x, c(1, 6, 11))
  expect_equal(out$d, c(0.1, 0.15, 0.2))
})
