context("Test forecasting context")

library(forecasting)
library(KernSmooth)
library(kde1d)
library(pracma)

test_that("kde lookup throws error", {
  expect_error(marginal_lookup(4), "Method selection*")
})

# -----------------------------------------------------------------------------
# Test containers for existing distribution estimation functions

test_that("probprecalc throws error", {
  expect_error(probprecalc(1:4), "Precalculated*")
})

test_that("probprecalc handles inputs with and without density precalculated", {
  fake_input <- structure(list("quantiles"=list(x=1:4, q=5:8)), class=c("fc_bma"))
  with_mock(calc_quantiles=function(y, ...) return(list(x=y$quantiles$x, q=y$quantiles$q, d=3:6)),
            expect_equal(probprecalc(fake_input), list(x=1:4, d=3:6, u=5:8)))
  fake_input <- structure(list("quantiles"=list(x=1:4, d=1:4, q=5:8)), class=c("fc_bma"))
  with_mock(calc_quantiles=function(y, ...) return(list(x=y$quantiles$x, q=y$quantiles$q, d=3:6)),
            expect_equal(probprecalc(fake_input), list(x=1:4, d=1:4, u=5:8)))
})

test_that("probempirical handles a given xmax", {
  with_mock(check_xmax=function(x, xmax) return(xmax),
    out <- probempirical(c(1,2, 4), xmax=5))
  expect_equal(out$x, c(0, 1, 2, 4, 5))
  expect_equal(out$u, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("probempirical replaces given xmax with actual data maximum", {
  with_mock(check_xmax=function(x, xmax) return(5.2),
            out <- probempirical(c(1, 5.2), xmax=5))
  expect_equal(out$x, c(0, 1,  5.2))
  expect_equal(out$u, c(0, 0.5, 1))
})

test_that("probempirical handles default max", {
  with_mock(check_xmax=function(x, xmax) return(5),
    out <- probempirical(c(1, 2, 4, 5)))
  expect_equal(out$x, c(0, 1,  2, 4, 5))
  expect_equal(out$u, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("probempirical ignores missing data", {
  with_mock(check_xmax=function(x, xmax) return(xmax),
    out <- probempirical(c(NaN, 1, 6), xmax=5))
  expect_equal(out$x, c(0, 1,  6))
  expect_equal(out$u, c(0, 0.5, 1))
  with_mock(check_xmax=function(x, xmax) return(5),
    out <- probempirical(c(1, 2, 4, NaN, 5)))
  expect_equal(out$x, c(0, 1,  2, 4, 5))
  expect_equal(out$u, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("probempirical handles both duplicate values. ", {
  with_mock(check_xmax=function(x, xmax) return(xmax),
    out <- probempirical(c(1, 2, 2, 4), xmax=5))
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
  expect_error(probtranskde(x, xmax="one"), "Bad input.*")
  expect_error(probtranskde(x, xmax=0), "Bad input.*")
})

test_that("transform_to_real_line throws errors", {
  # Negative value errors
  expect_error(transform_to_real_line(c(-1, 2, 4),  xmax=NaN, scale=0.5, zero_offset = 0.1), "Transformation*")
  expect_error(transform_to_real_line(c(-1, 2, 4),  xmax=10, scale=0.5, zero_offset = 0.1), "Transformation*")
  # zero errors
  expect_error(transform_to_real_line(c(0, 2, 4),  xmax=4, scale=0.9, zero_offset = 0), "Transformation*")
  expect_error(transform_to_real_line(c(0, 2, 4),  xmax=NaN, scale=0.9, zero_offset = 0), "Transformation*")
  # 1 errors
  expect_error(transform_to_real_line(c(0, 2, 4),  xmax=4, scale=1, zero_offset = 0.1), "Transformation*")
})

test_that("get_output_seq is correct.", {
  expect_equal(get_output_seq(x=c(0, 1, 3, 4), xmax=NaN, n.res=4, scaler=2), c(0.2, 0.4, 0.6, 0.8)*8)
  expect_equal(get_output_seq(x=c(1, 2, 3, 5),  xmax=3, n.res=4, scaler=2), c(0.2, 0.4, 0.6, 0.8))
})

test_that("get_output_seq ignores missing data.", {
  expect_equal(get_output_seq(x=c(0, 1, 3, 4, NA), xmax=NaN, n.res=4, scaler=2), c(0.2, 0.4, 0.6, 0.8)*8)
})

test_that("Scale_01 adjusts zeros for log transform.", {
  expect_equal(scale_01(x=c(0, 2, 4),  xmax=NaN, scale=0.5, zero_offset = 0.0001), c(0.0001, 2, 4))
})

test_that("Scale_01 rescales for probit transform.", {
  expect_equal(scale_01(x=c(0, 2, 4), xmax=10, scale=0.5, zero_offset = 0.0001), c(0.0001, 0.1, 0.2))
})

test_that("Scale_full rescales and adds 0 for probit transform.", {
  with_mock(cumtrapz=function(x, d) return(cumsum(d)),
            out <- scale_full(x01=c(0.1, 0.2, 0.3, 0.4, 0.5), d01=rep(1, times=5), xmax=10))
  expect_equal(out$x, seq(0, 10, by=2))
  expect_equal(out$u, seq(0, 1, by=0.2))
})

test_that("Scale_full adds 0 for log transform.", {
  with_mock(cumtrapz=function(x, d) return(cumsum(d)),
            out <- scale_full(x01=seq(0.1, 0.5, by=0.1), d01=rep(1, times=5), xmax=NaN))
  expect_equal(out$x, seq(0, 0.5, by=0.1))
  expect_equal(out$u, seq(0, 1, by=0.2))
})

test_that("Scale_full doesn't add 0 when it's already there.", {
  with_mock(cumtrapz=function(x, d) return(cumsum(d)),
            out <- scale_full(x01=seq(0, 0.5, by=0.1), d01=rep(1, times=6), xmax=NaN))
  expect_equal(out$x, seq(0, 0.5, by=0.1))
})

test_that("check_xmax is correct", {
  expect_equal(check_xmax(c(0, 1, 2, 3), xmax=4), 4)
  expect_equal(check_xmax(c(0, 1, 2, 3), xmax=2), 3)
})

test_that("check_xmax ignores NA's", {
  expect_equal(check_xmax(c(0, 1, 2, 3, NaN), xmax=2), 3)
})
