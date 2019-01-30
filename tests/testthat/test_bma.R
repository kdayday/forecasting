context("Test forecasting context")

library(forecasting)

test_that("Probability of clipping calculation is correct.", {
  expect_equal(get_poc(FCST=NA, A0=0.1, A1=1.2, A2=1), NA)
  expect_equal(get_poc(FCST=0.5, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-0.7))) # 0.1+1.2*0.5
  expect_equal(get_poc(FCST=1, A0=0.1, A1=1.2, A2=1), 1/(1+exp(-2.3))) # 0.1+1.2*1+1
})

test_that("get_z is correct.", {
  expect_equal(get_z(OBS=1, PoC=0.4, db=0.01, w=0.5), 0.2) # 0.5 * 0.4 = 0.2
  expect_equal(get_z(OBS=0.8, PoC=0.4, db=0.01, w=0.5), 0.003) # 0.5 * 0.6 * 0.01 =0.003
})

test_that("get_z handles missing members", {
  expect_true(is.na(get_z(OBS=NA, PoC=0.4, db=0.01, w=0.5)))
  expect_true(is.na(get_z(OBS=0.8, PoC=NA, db=0.01, w=0.5)))
})

test_that("get_beta_density is correct", {
  # mu = 0.1 + 0.4=0.5
  # sigma =0.1 + 0.1 = 0.2
  # gamma = 5.25
  with_mock(dbeta=function(...) return(sum(...)),
            expect_equal(get_beta_density(OBS=0.75, FCST=0.5, B0=0.1, B1=0.8, C0=0.1, C1=0.2), 6))
})

test_that("get_beta_density throws error", {
  # Get negative gamma calculation -- PLACEHOLDER UNTIL I DOUBLE-CHECK/FIX THIS
  expect_error(get_beta_density(OBS=0.75, FCST=0.5, B0=0.05, B1=0.1, C0=0.5, C1=0.8), "Non-positive*")
})

test_that("get_beta_density is returns NA's", {
  expect_equal(get_beta_density(OBS=NA, FCST=0.5, B0=0.1, B1=0.8, C0=0.1, C1=0.2), NA)
  expect_equal(get_beta_density(OBS=0.5, FCST=NA, B0=0.1, B1=0.8, C0=0.1, C1=0.2), NA)
  expect_equal(get_beta_density(OBS=1, FCST=0.5, B0=0.1, B1=0.8, C0=0.1, C1=0.2), NA)
})

test_that("e_step array handling is correct.", {
  FCST <- array(data=2*(1:8), dim = c(2,2,2))
  B0 <- array(data=-1*(1:8), dim = c(2,2,2))
  B1 <- array(data=1:8, dim = c(2,2,2))
  C0 <- 1
  C1 <- 2
  w <- c(10, 20)
  PoC <- array(data=3*(1:8), dim = c(2,2,2))
  OBS <- matrix(1:4, ncol=2)

  with_mock(get_beta_density=function(OBS, FCST, B0, B1, C0, C1) {return(sum(OBS, FCST, B0, B1, C0, C1))}, # 1+ 2-1+1+(1+2) = 6; c(6, 9, 12, 15, 14, 17, 20, 23)
            get_z= function(OBS, PoC, db, w) return(db-PoC+OBS+w), # c(3, 3, 3, 3, -1, -1, -1, -1) -> c(4, 5, 6, 7, 0, 1, 2, 3) -> c(14, 15, 16, 17, 20, 21, 22, 23)
            OUT <- e_step(w, C0, C1, OBS, FCST, B0, B1, PoC))
  expect_equal(OUT$z, array(c(14/34, 15/36, 16/38, 17/40, 20/34, 21/36, 22/38, 23/40), dim=c(2,2,2)))
  expect_equal(OUT$sumz, matrix(c(34, 36, 38, 40), ncol=2))
})

test_that("em_subfunction is correct.", {
  with_mock(e_step=function(...) {return(list(z=(array(c(0.1, 0.2, 0.3, NA, seq(0.2, 0.8, by = 0.2), rep(NA, 4)), dim=c(2, 2, 3)))))},
            optim=function(...) {return(list(convergence=0, par=c(2, 3)))},
            OUT <- em_subfunction(FCST=array(1:12, dim = c(2,2,3)), OBS=NA, PoC=NA, B0=NA, B1=NA, C0=4, C1=6, w=c(1, 0.7, 0)))
  expect_equal(OUT$w, c(0.2, 0.5, 0))
  expect_equal(OUT$C0, 2)
  expect_equal(OUT$C1, 3)
  expect_equal(abs(OUT$error), c(0.8, 0.2, 0, 2, 3))
})

