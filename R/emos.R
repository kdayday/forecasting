#' Fit EMOS coefficients a, b_k, c, d for truncated normal
#' b_k, c, and d are parameterized as squares to control non-negativity of regression and scale parameters, like the "square" coefRule and varRule in the ensembleMOS package.
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
    par_init <- list(a=start_coefs[1], b=start_coefs[-1], c=5, d=1)
  }

  # Square-root b, c, d parameters
  par <- c(sqrt(par_init$c), sqrt(par_init$d), par_init$a, sqrt(abs(par_init$b)))

  lb <- c(0, 0, rep(-Inf, times=length(par)-2)) # Force variance parameters to be lower bounded at 0
  opt <- optim(par, tnorm_crps_obj, obs=tel, ens=ens, max_power=max_power, method="L-BFGS-B", lower=lb, control=list(maxit=1e7))

  # Re-square b, c, d parameters
  return(list(a=opt$par[3], b=opt$par[-(1:3)]^2, c=opt$par[1]^2, d=opt$par[2]^2))
}


#' Subfunction to get average crps over the training period
#' @param par a vector of parameters to optimize over: c, d, a, b1, ... bn
#' @param obs a vector of observations of length m
#' @param ens a [time x member] (i.e., [m x n]) matrix of the ensemble forecasts for each observation
#' @param max_power Site AC power, for upper limit of truncated normal
tnorm_crps_obj <- function(par, obs, ens, max_power) {
  # Re-square parameters
  C <- par[1]^2
  D <- par[2]^2
  A <- par[3]
  B <- par[-(1:3)]^2

  if (length(B) != ncol(ens)) stop("Mismatch between number of parameters and number of ensemble members.")

  mn <- apply(X=ens, MARGIN = 1, FUN=function(e) return(A + sum(B*e)))
  S2 <- apply(X=ens, MARGIN = 1, FUN=var)
  std_dev <- C + D*S2
  crps <- mean(scoringRules::crps_tnorm(obs, location = mn, scale = std_dev, lower=0, upper=max_power), na.rm=T)

  return(crps)
}



