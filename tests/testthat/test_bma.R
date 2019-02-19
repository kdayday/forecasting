context("Test forecasting context")

library(forecasting)

test_that("Probability of clipping calculation is correct.", {
  expect_equal(get_poc(FCST=NA, A0=0.1, A1=1.2, A2=1), NA)
  expect_equal(get_poc(FCST=0.5, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-0.7))) # 0.1+1.2*0.5
  expect_equal(get_poc(FCST=1, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-2.3))) # 0.1+1.2*1+1
})

test_that("get_z is correct, including clipping tolerance", {
  expect_equal(get_z(OBS=1-1e-7, PoC=0.4, db=0.01, w=0.5, tol=1e-6), 0.2) # 0.5 * 0.4 = 0.2
  expect_equal(get_z(OBS=1-1e-5, PoC=0.4, db=0.01, w=0.5, tol=1e-6), 0.003) # 0.5 * 0.6 * 0.01 =0.003
})

test_that("get_z handles missing members", {
  expect_true(is.na(get_z(OBS=NA, PoC=0.4, db=0.01, w=0.5, tol=1e-6)))
  expect_true(is.na(get_z(OBS=0.8, PoC=NA, db=0.01, w=0.5, tol=1e-6)))
})

test_that("e_step array handling is correct.", {
  FCST <- array(data=2*(1:8), dim = c(2,2,2))
  B0 <- array(data=-1*(1:8), dim = c(2,2,2))
  B1 <- array(data=1:8, dim = c(2,2,2))
  C0 <- 1
  w <- c(10, 20)
  PoC <- array(data=3*(1:8), dim = c(2,2,2))
  OBS <- matrix(1:4, ncol=2)

  with_mock(get_rho= function(FCST, B0, B1, ...) return(sum(FCST, B0, B1)),
            get_gamma = function(rho, C0) return(sum(rho, C0)), # 3:2:17
            dbeta_gamma_rho = function(OBS, gammas, rhos) return(sum(OBS, gammas, rhos)), # 6, 11, 16, 21, 22, 27, 32, 37
            get_z= function(OBS, PoC, db, w, tol) return(db-PoC+OBS+w), #  c(4, 7, 10, 13, 8, 11, 14, 17 ) -> c(14, 17, 20, 23, 28, 31, 34, 37)
            OUT <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform=NA, tol=NA))
  expect_equal(OUT$z, array(c(14/42, 17/48, 20/54, 23/60, 28/42, 31/48, 34/54, 37/60), dim=c(2,2,2)))
  expect_equal(OUT$sumz, matrix(c(42, 48, 54, 60), ncol=2))
})

test_that("e_step sumz handles NaN's.", {
  OBS <-array(c(NA, 1, 1, NA, NA, 2), dim=c(3, 1, 2))
  with_mock(get_rho= function(...) return(NA), get_gamma = function(...) return(NA),
            dbeta_gamma_rho = function(...) return(NA),
            get_z= function(x, ...) return(x),
            OUT <- e_step(w=NA, C0=NA, OBS=OBS, FCST=OBS, B0=NA, B1=NA, PoC=NA, B_transform=NA, tol=NA))
  expect_equal(OUT$z, array(c(NA, 1, 1/3, NA, NA, 2/3), dim=c(3, 1, 2)))
  expect_equal(OUT$sumz, matrix(c(NA, 1, 3), ncol=1))
})

test_that("em_subfunction is correct.", {
  with_mock(e_step=function(...) {return(list(z=(array(c(0.1, 0.2, 0.3, NA, seq(0.2, 0.8, by = 0.2), rep(NA, 4)), dim=c(2, 2, 3)))))},
            optim=function(...) {return(list(convergence=0, par=2))},
            OUT <- em_subfunction(FCST=array(1:12, dim = c(2,2,3)), OBS=NA, PoC=NA, B0=NA, B1=NA, C0=4, w=c(1, 0.7, 0), tol=NA))
  expect_equal(OUT$w, c(0.2, 0.5, 0))
  expect_equal(OUT$C0, 2)
  expect_equal(abs(OUT$error), c(0.8, 0.2, 0, 2))
})

test_that("beta1_ens_modesl throws errors", {
  expect_error(beta1_ens_models(tel=c(1, 2, NA), ens=matrix((1:9)/9, ncol=3)), "Telemetry*")
  expect_error(beta1_ens_models(tel=c(1, 2, NA)/3, ens=matrix((1:9), ncol=3)), "All forecasts*")
  expect_error(beta1_ens_models(tel=c(1, 2, NA)/3, ens=matrix((1:9)/9, ncol=1)), "Must*")
})

test_that("beta1_ens_modesl coefficient data look-up is correct", {
  fake_lr <- function(e, tel, ...) {
    coeffs <- data.frame(Estimate=c(sum(e), sum(tel)), P=c(e[1],e[3]), row.names = c("(Intercept)", "x"))
    names(coeffs) <- c("Estimate", "Pr(>|z|)") # Overwrite with special character names
    return(list(coefficients=coeffs, aic=10))
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

test_that("get_lm clipping tolerance is correct", {
  with_mock(lm= function(form, data, ...) return(data),
            summary = function(x) return(x),
    OUT <- get_lm(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1.001, 0.999), form=NA, B_transform=function(x) return(x), tol=1e-2)
  )
  expect_equal(dim(OUT), c(1,2))
})

test_that("get_lr clipping tolerance is correct", {
  with_mock(glm= function(form, data, ...) return(data),
            summary = function(x) return(x),
            OUT <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1.001, 0.999), form=NA, A_transform=function(x) return(x), tol=1e-2)
  )
  expect_equal(OUT[,'y'], c(0, 1, 1))
})

test_that("get_lm handles transform", {
  with_mock(lm= function(form, data, ...) return(data),
            summary = function(x) return(x),
            OUT_transform <- get_lm(fc=c(0.2), tel=c(0.5), form=NA, B_transform=function(x) return(x+0.1), tol=1e-2),
            OUT_no_transform <- get_lm(fc=c(0.2), tel=c(0.5), form=NA, B_transform=NA, tol=1e-2)
  )
  expect_equal(OUT_transform[,'x'], c(0.3))
  expect_equal(OUT_no_transform[,'x'], c(0.2))
})

test_that("get_lr handles transform", {
  with_mock(glm= function(form, data, ...) return(data),
            summary = function(x) return(x),
            OUT_transform <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1.001, 0.999), form=NA, A_transform=function(x) return(x+0.1), tol=1e-2),
            OUT_no_transform <- get_lr(fc=c(0.2, 0.4, 0.6), tel=c(0.5, 1.001, 0.999), form=NA, A_transform=NA, tol=1e-2)
  )
  expect_equal(OUT_transform[,'x'], c(0.3, 0.5, 0.7))
  expect_equal(OUT_no_transform[,'x'], c(0.2, 0.4, 0.6))
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


test_that("get_gamma handles NA's", {
  expect_true(is.na(get_gamma(mu=NA, C0=0.1)))
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
            expect_equal(get_log_lik(C0=NA, w=NA, OBS=NA, FCST=NA, B0=NA, B1=NA, PoC=NA, B_transform=NA, tol=NA), 7)
  )
})
