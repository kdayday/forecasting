#' Fit EMOS coefficients a, b_k, c, d for truncated normal
#' @param tel Vector of training telemetry data
#' @param ens Matrix of training ensemble member data [time x member]
#' @param par_init Initial parameter values as a list with a, b, c, d
#' @param max_power site AC rating, for upper limit of truncated normal
#'
emos_model <- function(tel, ens, max_power, par_init=NA) {
  if (!(length(tel) == nrow(ens))) stop("Must have same number of telemetry and forecast time-points.")

  # Some defaults are referenced to fitMOStruncnorm function in ensembleMOS
  if (all(is.na(par_init))) {
    # Start off with a simple linear model
    start_coefs <- unname(lm(tel~ens)$coef)
    # C, D, A, B
    par <- c(5, 1, start_coefs)
  } else {par <- c(par_init$c, par_init$d, par_init$a, par_init$b)}

  lb <- c(0, 0, rep(-Inf, times=length(par)-2)) # Force variance parameters to be lower bounded at 0
  opt <- optim(par, tnorm_crps_obj, obs=tel, ens=ens, max_power=max_power, method="L-BFGS-B", lower=lb, control=list(maxit=1e7))

  return(list(a=opt$par[3], b=opt$par[-(1:3)], c=opt$par[1], d=opt$par[2]))
}


#' Subfunction to get average crps over the training period
#' @param par a vector of parameters to optimize over: c, d, a, b1, ... bn
#' @param obs a vector of observations of length m
#' @param ens a [time x member] (i.e., [m x n]) matrix of the ensemble forecasts for each observation
#' @param max_power Site AC power, for upper limit of truncated normal
tnorm_crps_obj <- function(par, obs, ens, max_power) {
  C <- par[1]
  D <- par[2]
  A <- par[3]
  B <- par[-(1:3)]

  if (length(B) != ncol(ens)) stop("Mismatch between number of parameters and number of ensemble members.")

  mn <- apply(X=ens, MARGIN = 1, FUN=function(e) return(A + sum(B*e)))
  S2 <- apply(X=ens, MARGIN = 1, FUN=var)
  std_dev <- C + D*S2
  crps <- mean(scoringRules::crps_tnorm(obs, location = mn, scale = std_dev, lower=0, upper=max_power), na.rm=T)

  if (!is.finite(crps)) browser()
  return(crps)
}



