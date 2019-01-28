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
