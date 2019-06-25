context("Test forecasting context")

library(forecasting)
library(logistf)

test_that("Probability of clipping calculation is correct.", {
  expect_equal(get_poc(FCST=NA, A0=0.1, A1=1.2, A2=1), NA)
  expect_equal(get_poc(FCST=0.5, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-0.7))) # 0.1+1.2*0.5
  expect_equal(get_poc(FCST=1, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-2.3))) # 0.1+1.2*1+1
})

test_that("get_weighted_probability is correct, including clipping threshold", {
  # Clipped
  with_mock(dbeta = function(obs, a, b) return(a + b),
            expect_equal(get_weighted_probability(OBS=0.99, PoC=0.4, gamma=0.1, rho=0.4, w=0.5, percent_clipping_threshold=0.99), 20)) # 0.5 * 0.4/0.01 = 20
  # Unclipped
  with_mock(dbeta =function(obs, a, b) return(a+b), # db = 0.01
            pbeta =function(obs, a, b) return(0.8),
    expect_equal(get_weighted_probability(OBS=0.989, PoC=0.4, gamma=0.01, rho=0.4, w=0.5, percent_clipping_threshold=0.99), 0.00375)) # 0.5 * 0.6 * 0.01/0.8 =0.00375
})

test_that("get_weighted_probability handles missing members", {
  expect_true(is.na(get_weighted_probability(OBS=NA, PoC=0.4, gamma=0.01, rho=0.4, w=0.5, percent_clipping_threshold=0.99)))
  expect_true(is.na(get_weighted_probability(OBS=0.8, PoC=NA, gamma=0.01, rho=0.4, w=0.5, percent_clipping_threshold=0.99)))
})

test_that("get_weighted_probability returns missing for observataions at 0", {
  expect_true(is.na(get_weighted_probability(OBS=0, PoC=0.4, gamma=0.01, rho=0.4, w=0.5, percent_clipping_threshold=0.99)))
})

test_that("e_step array handling is correct.", {
  FCST <- array(data=2*(1:8), dim = c(2,2,2))
  B0 <- array(data=-1*(1:8), dim = c(2,2,2))
  B1 <- array(data=1:8, dim = c(2,2,2))
  C0 <- 1
  w <- c(10, 20)
  PoC <- array(data=3*(1:8), dim = c(2,2,2))
  OBS <- array(rep(1:4, times=2), dim=c(2,2,2))

  with_mock(get_rho= function(FCST, B0, B1, ...) return(sum(FCST, B0, B1)),
            get_gamma = function(rho, C0) return(sum(rho, C0)), # 3:2:17
            # 6, 11, 16, 21, 22, 27, 32, 37 ->  c(4, 7, 10, 13, 8, 11, 14, 17 ) -> c(14, 17, 20, 23, 28, 31, 34, 37)
            get_weighted_probability= function(OBS, PoC, gamma, rho, w, percent_clipping_threshold) return(sum(OBS, gamma, rho)-PoC+OBS+w),
            OUT <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform=NA, percent_clipping_threshold=NA))
  expect_equal(OUT$z, array(c(14/42, 17/48, 20/54, 23/60, 28/42, 31/48, 34/54, 37/60), dim=c(2,2,2)))
  expect_equal(OUT$sumz, matrix(c(42, 48, 54, 60), ncol=2))
})

test_that("e_step sumz handles NaN's.", {
  OBS <-array(c(NA, 1, 1, NA, NA, 2), dim=c(3, 1, 2))
  with_mock(get_rho= function(...) return(NA), get_gamma = function(...) return(NA),
            get_weighted_probability= function(x, ...) return(x),
            OUT <- e_step(w=NA, C0=NA, OBS=OBS, FCST=OBS, B0=NA, B1=NA, PoC=NA, B_transform=NA, percent_clipping_threshold=NA))
  expect_equal(OUT$z, array(c(NA, 1, 1/3, NA, NA, 2/3), dim=c(3, 1, 2)))
  expect_equal(OUT$sumz, matrix(c(NA, 1, 3), ncol=1))
})

test_that("em_subfunction is correct.", {
  with_mock(e_step=function(...) {return(list(z=(array(c(0.1, 0.2, 0.3, NA, seq(0.3, 0.9, by = 0.2), rep(NA, 4)), dim=c(2, 2, 3)))))},
            optimize=function(...) {return(list(maximum=2))},
            OUT <- em_subfunction(FCST=array(1:12, dim = c(2,2,3)), OBS=NA, PoC=NA, B0=NA, B1=NA, C0=4, w=c(0.2, 0.8, 0),
                                  percent_clipping_threshol=NA, count=1, CM2.iter=1))
  expect_equal(OUT$w, c(0.25, 0.75, 0))
  expect_equal(sum(OUT$w, na.rm=T), 1)
  expect_equal(OUT$C0, 2)
  expect_equal(abs(OUT$error), c(0.05, 0.05, 0, 2))
})

test_that("beta1_ens_modesl throws errors", {
  expect_error(beta1_ens_models(tel=c(1, 1.01, NA), ens=matrix((1:9)/9, ncol=3)), "Telemetry*")
  expect_error(beta1_ens_models(tel=c(1, 2, NA)/3, ens=matrix((1:9), ncol=3)), "All forecasts*")
  expect_error(beta1_ens_models(tel=c(1, 2, NA)/3, ens=matrix((1:9)/9, ncol=1)), "Must*") # telemetry vector wrong length
  # Does not check for correct order of dimensions
  expect_error(beta1_ens_models(tel=matrix((1:6)/9, ncol=3), ens=matrix((1:9)/9, ncol=3)), "Must*") # telemetry matrix wrong size
  expect_warning(beta1_ens_models(tel=c(1, 2, NA)/3, ens=matrix((1:9)/9, ncol=3), percent_clipping_threshold=1.1), "Percent*")
})

test_that("beta1_ens_modesl coefficient data look-up is correct", {
  fake_lr <- function(e, tel, ...) {
    coeffs <- c(sum(e), sum(tel))
    names(coeffs) <- c("(Intercept)", "x")
    probs <- c(e[1],e[3])
    names(probs) <- c("(Intercept)", "x")
    return(list(coefficients=coeffs, prob=probs, aic=10))
  }
  fake_lm <- function(e, tel, ...) {
    coeffs <- data.frame(Estimate=c(sum(e), sum(tel)), P=c(2,3), row.names = c("(Intercept)", "x"))
    names(coeffs) <- c("Estimate", "Pr(>|t|)") # Overwrite with special character names
    return(list(coefficients=coeffs, r.squared=0.5))
  }
  with_mock(get_lr= fake_lr, get_lm=fake_lm,
            em = function(...) return(list(C0=NA, C1=NA, w=NA, log_lik=NA)),
            OUT <- beta1_ens_models(tel=c(0.25, 0.5, 0.75), ens=matrix(seq(0.1, 0.9, by=0.1), ncol=3))
  )
  expect_equal(OUT$A0, c(0.6, 1.5, 2.4))
  expect_equal(OUT$A1, c(1.5, 1.5, 1.5))
  expect_equal(OUT$fit_statistics, data.frame("A0 p-value"=c(0.1, 0.4, 0.7),
                                              "A1 p-value"=c(0.3, 0.6, 0.9),
                                              "A AIC"=c(10,10,10),
                                              "B0 p-value"=c(2, 2, 2),
                                              "B1 p-value"=c(3, 3, 3),
                                              "B r-squared"=c(0.5, 0.5, 0.5)))
})


test_that("beta1_ens_modesl coefficient data look-up handles missing fields", {
  fake_lr <- function(e, tel, ...) {
    coeffs <- c(sum(e), sum(tel))
    names(coeffs) <- c("(Intercept)", "x")
    probs <- c(e[1],e[3])
    names(probs) <- c("(Intercept)", "x")
    return(list(coefficients=coeffs, prob=probs, aic=10))
  }
  fake_lm <- function(e, tel, ...) {
    coeffs <- data.frame(Estimate=numeric(), p=numeric())
    names(coeffs) <- c("Estimate", "Pr(>|t|)") # Overwrite with special character names
    return(list(coefficients=coeffs, r.squared=0))
  }
  with_mock(get_lr= fake_lr, get_lm=fake_lm,
            em = function(...) return(list(C0=NA, C1=NA, w=NA, log_lik=NA)),
            OUT <- beta1_ens_models(tel=c(0.25, 0.5, 0.75), ens=matrix(seq(0.1, 0.9, by=0.1), ncol=3))
  )
  expect_equal(OUT$B0, c(0, 0, 0))
  expect_equal(OUT$B1, c(NA, NA, NA))
  expect_equal(OUT$fit_statistics, data.frame("A0 p-value"=c(0.1, 0.4, 0.7),
                                              "A1 p-value"=c(0.3, 0.6, 0.9),
                                              "A AIC"=c(10,10,10),
                                              "B0 p-value"=c(NA, NA, NA),
                                              "B1 p-value"=c(NA, NA, NA),
                                              "B r-squared"=c(0, 0, 0)))
})

test_that("get_lm clipping threshold is correct", {
  with_mock(lm= function(form, data, ...) return(data),
            summary = function(x) return(x),
    OUT <- get_lm(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1, 0.99), form=NA, B_transform=function(x) return(x), percent_clipping_threshold=0.99)
  )
  expect_equal(dim(OUT), c(1,2))
})

test_that("get_lm handles missing or no unclipped data", {
  OUT <- get_lm(fc=c(0.2, 0.4, 0.6), tel=c(NA, 0.95, 0.99), form=NA, B_transform=function(x) return(x), percent_clipping_threshold=0.95)
  expect_true(all(is.na(OUT$coefficients)))
})

test_that("get_lr clipping threshold is correct", {
  with_mock(logistf= function(form, data, ...) return(data),
            extractAIC = function(...) return(NA),
            summary = function(x) return(x),
            OUT <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1, 0.99), form=NA, A_transform=function(x) return(x), percent_clipping_threshold=0.99)
  )
  expect_equal(OUT[['y']], c(0, 1, 1))
})

test_that("get_lr handles complete single outcome value", {
  OUT <- get_lr(fc=c(0.2, 0.4, 0.6, NA, 0), tel=c(0, 0, NA, 0.1, 0.1), form=NA, A_transform=function(x) return(x), percent_clipping_threshold=0.99)
  expect_equal(unname(OUT$coefficients), c(-Inf, 0))
  expect_true(is.na(OUT$aic))

  OUT <- get_lr(fc=c(0.2, 0.4, 0.6, NA, 0), tel=c(1, 0.999, NA, 0.9, 0.9), form=NA, A_transform=function(x) return(x), percent_clipping_threshold=0.99)
  expect_equal(unname(OUT$coefficients), c(+Inf, 0))
})

test_that("get_lm handles transform", {
  with_mock(lm= function(form, data, ...) return(data),
            summary = function(x) return(x),
            OUT_transform <- get_lm(fc=c(0.2), tel=c(0.5), form=NA, B_transform=function(x) return(x+0.1), percent_clipping_threshold=0.99),
            OUT_no_transform <- get_lm(fc=c(0.2), tel=c(0.5), form=NA, B_transform=NA, percent_clipping_threshold=0.99)
  )
  expect_equal(OUT_transform[,'x'], c(0.3))
  expect_equal(OUT_no_transform[,'x'], c(0.2))
})

test_that("get_lr handles transform", {
  with_mock(logistf= function(form, data, ...) return(data),
            extractAIC = function(...) return(NA),
            summary = function(x) return(x),
            OUT_transform <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1, 0.99), form=NA, A_transform=function(x) return(x+0.1), percent_clipping_threshold=0.99),
            OUT_no_transform <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1, 0.99), form=NA, A_transform=NA, percent_clipping_threshold=0.99)
  )
  expect_equal(OUT_transform[['x']], c(0.3, 0.5, 0.7))
  expect_equal(OUT_no_transform[['x']], c(0.2, 0.4, 0.6))
})

test_that("get_rho handles transforms correctly", {
  expect_equal(get_rho(FCST=0.5, B0=0.25, B1=1, B_transform=NA), 0.75)
  expect_equal(get_rho(FCST=0.5, B0=0.25, B1=1, B_transform=function(x) return(x-0.25)), 0.5)
})

test_that("get_rho truncates a touch below 1",{
  expect_equal(get_rho(FCST=2, B0=1, B1=0.5, B_transform=NA), 1-1e-6)
})

test_that("get_rho handles NA's", {
  expect_true(is.na(get_rho(FCST=NA, B0=1, B1=0.5, B_transform=NA)))
  expect_true(is.na(get_rho(FCST=0.5, B0=NA, B1=NA, B_transform=NA)))
})

test_that("get_gamma calculation is correct", {
  # Example under the theoretical limit
  expect_equal(get_gamma(mu=0.7, C0=0.1), 0.21/0.084 - 1) # sigma = sqrt(0.1-0.4*0.04), gamma = 0.21/0.984 - 1
  # Example over the theoretical limit
  expect_equal(get_gamma(mu=0.7, C0=0.25), 0.21/(sqrt(0.21)-1e-6)^2 - 1) # sigma = sqrt(0.21)-1e-6
})


test_that("get_gamma handles NA's and mean on boundary", {
  expect_true(is.na(get_gamma(mu=NA, C0=0.1)))
  expect_true(is.na(get_gamma(mu=0, C0=0.1)))
})

test_that("dbeta_gamma_rho is correct", {
  expect_equal(dbeta_gamma_rho(x=0.5, g=2+2, rho=2/(2+2)), 1.5) # for alpha=beta=2
})

test_that("dbeta_gamma_rho handles NA's", {
  expect_true(is.na(dbeta_gamma_rho(x = NA, g=4, rho = 0.5)))
  expect_true(is.na(dbeta_gamma_rho(x = 0.5, g=NA, rho = NA)))
})

test_that("get_log_lik is correct",{
  with_mock(e_step=function(...) return(list(sumz=c(NA, exp(2), exp(5)))),
            expect_equal(get_log_lik(C0=NA, w=NA, OBS=NA, FCST=NA, B0=NA, B1=NA, PoC=NA, B_transform=NA, percent_clipping_threshold=NA), 7)
  )
})
