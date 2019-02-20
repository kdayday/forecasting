

# Fit BMA coefficients for beta distributions with discrete component at 1 for each ensemble member
# Data should be pre-processed and normalized to [0,1] so that clipped values are exactly 1
#' @param tel Vector of training telemetry data on [0,1]
#' @param ens Matrix of training ensemble member data [time x member] on [0,1]
#' @param lr_formula Formula in terms of x,y for logistic regression model, defaults to "y ~ x". Requires a negative x intercept to model PoC < 0.5.
#' @param A_transform A function for transforming forecast data before logistic regression to get a's (optional)
#' @param lm_formula Formula in terms of x,y for linear regression model, defaults to "y ~ x + 0"
#' @param B_transform A function for transforming forecast data before linear regression to get b's (optional)
#' @param ... Optional arguments to pass to EM algorithm
#' @return A formula for a discrete-continuous mixture model with beta distribution
beta1_ens_models <- function(tel, ens, lr_formula= y ~ x, A_transform=NA, lm_formula= y ~ x + 0, B_transform=NA, tol.clip=1e-6, ...) {
  if (any(tel < 0, na.rm=T) | any(tel-1>tol.clip, na.rm=T)) stop('Telemetry must be normalized to [0,1] to apply beta model.')
  if (any(ens < 0, na.rm=T) | any(ens-1>tol.clip, na.rm=T)) stop('All forecasts must be normalized to [0,1] to apply beta model.')
  if (length(tel) != dim(ens)[1]) stop("Must have same number of telemetry and forecast time-points.")

  # 1. Logistic regression for a's
  # Returns a list of model summaries, one per member
  mem_discrete_models <-lapply(seq_len(dim(ens)[2]), function(i) get_lr(ens[,i], tel=tel, form=lr_formula, A_transform = A_transform, tol.clip=tol.clip))
  A0 <- sapply(mem_discrete_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Estimate"], 0)))
  A1 <- sapply(mem_discrete_models, function(m) return(m$coefficients["x", "Estimate"]))

  # 2. linear regression for b's.
  mem_mean_models <- lapply(seq_len(dim(ens)[2]), function(i) get_lm(ens[,i], tel=tel, form=lm_formula, B_transform = B_transform, tol.clip=tol.clip))
  B0 <- sapply(mem_mean_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Estimate"], 0)))
  B1 <- sapply(mem_mean_models, function(m) return(m$coefficients["x", "Estimate"]))

  fit_statistics <- data.frame("A0 p-value"=sapply(mem_discrete_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Pr(>|z|)"], NA))),
                               "A1 p-value"=sapply(mem_discrete_models, function(m) return(m$coefficients["x", "Pr(>|z|)"])),
                               "A AIC"=sapply(mem_discrete_models, function(m) return(m$aic)),
                               "B0 p-value"=sapply(mem_mean_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Pr(>|t|)"], NA))),
                               "B1 p-value"=sapply(mem_mean_models, function(m) return(m$coefficients["x", "Pr(>|t|)"])),
                               "B r-squared"=sapply(mem_mean_models, function(m) return(m$r.squared)))

  # 3. ME algorithm for w's and c's
  # In future, the above could be expanded to have unique site values. For now, matrices are expanded to arrays assuming a "single" site (i.e., global values).
  ntime <- dim(ens)[1]
  nens <- dim(ens)[2]
  array_dims <- c(ntime, 1, nens)
  tmp <- em(FCST=array(ens, dim=array_dims), OBS=array(tel, dim=c(ntime, 1)), A0=array(rep(A0, each=ntime), dim=array_dims),
            A1=array(rep(A1, each=ntime), dim=array_dims), A2=0, B0=array(rep(B0, each=ntime), dim=array_dims),
            B1=array(rep(B1, each=ntime), dim=array_dims), A_transform=A_transform, B_transform=B_transform, tol.clip=tol.clip, ...)

  return(list(A0=A0, A1=A1, B0=B0, B1=B1, C0=tmp$C0, w=tmp$w, fit_statistics=fit_statistics, log_lik=tmp$log_lik,
              A_transform=A_transform, B_transform=B_transform, em_count=tmp$count, em_error=tmp$max_error))
}


# Do logistic regression to get a single ensemble member's 'a' coefficients
# Telemetry within given tolerance of 1 are assumed to be clipped
# glm and lm omit NA values by default
# KEEP AN EYE ON SEPARATION PROBLEMS, AND JUST GET RID OF DISCRETE COMPONENT BASED ON ALL(CLIPPED==0, NA.RM=T) IF NEED BE.
#' @param fc Vector of training forecast data on [0,1]
#' @param tel Vector of telemetry data on [0,1]
#' @param form Formula in terms of x,y for logistic regression model
#' @param A_transform A function for pre-transforming the forecast data (optional), e.g. function(X) return(log(x)+1)
#' @param tol.clip Tolerance for determining if clipping is occurring
#' @return Summary of the glm model
get_lr <- function(fc, tel, form, A_transform, tol.clip){
  clipped <- as.integer(abs(tel-1) < tol.clip)
  if (typeof(A_transform)=="closure") fc <- sapply(fc, FUN = A_transform)
  return(summary(glm(form, family=binomial(link='logit'), data=data.frame(x=fc, y=clipped))))
}

# Get the linear model for a single ensemble member's 'b' coefficients
# Telemetry within given tolerance of 1 are assumed to be clipped
# glm and lm omit NA values by default
#' @param fc Vector of training forecast data on [0,1]
#' @param tel Vector of telemetry data on [0,1]
#' @param form Formula in terms of x,y for logistic regression model
#' @param B_transform A function for pre-transforming the forecast data (optional), e.g. function(X) return(log(x)+1)
#' @param tol.clip Tolerance for determining if clipping is occurring
#' @return Summary of the lm model
get_lm <- function(fc, tel, form, B_transform, tol.clip){
  unclipped <- abs(tel-1) > tol.clip
  if (typeof(B_transform)=="closure") fc <- sapply(fc, FUN = B_transform)
  return(summary(lm(form, data=data.frame(x=fc[unclipped], y=tel[unclipped]))))
}

# Expectation-maximization function, modified from code courtesy of Will Kleiber
em <- function(FCST, OBS, A0, A1, A2, B0, B1, A_transform, B_transform, tol.clip, C0=0.06, eps=1e-005, maxiter=1000, CM2.iter=50, start.w=NULL)

    # MODEL for one forecast : y_s is solar power, logit P(y_s = 1 | f) = a0 + a1 f
    #                          solar power level, conditional on it being less than rated power (i.e., 1) is beta distributed
    #                          with mean = b0 + b1 f, standard deviation = -C0/0.25*(mean-0.5)^2 + C0

    # Inputs:
    #  FCST        array of dimension (length of training period)x(number of stations)x(number of ensemble members)
    #              where each matrix (FCST[,,k]) has a ensemble member's forecast for these sites and times
#  OBS         matrix of dimension (length of training period)x(number of stations)
#              of observations at these sites and times
#  XN          for X=A,B and N=0,1,2, array of dimension (length of training period)x(number of stations)x
#              (number of ensemble members) of estimates of XN (usually posterior mean for A0,1,2 and
#              least squares from regression for B0,1)
#               Can set A2's to 0's to ignore the indicator function aspect

#  C0 starting estimate of C0, which is assumed constant across sites and members (equal variances among ensemble member)
#  eps     stopping criterion
#  maxiter maximum number of EM iterations allowed
#  CM2.iter CM-2 step of EM algorithm will be run every CM2.iter iterations of the M and CM-1 steps
#  start.w initial values for the weights (optional)

{
  # set intial weights
  w <- get_initial_weights(start.w, avail=apply(FCST, MARGIN = 3, FUN = function(x) {!all(is.na(x))}))

  # initialize values to get into loop
  error <- 1
  count <- 0

  # Precalculate probability of clipping and beta density
  PoC <- array(mapply(get_poc, FCST, A0, A1, A2, MoreArgs=list(A_transform=A_transform)), dim(FCST))

  # main EM algorithm
  while((max(abs(error)) > eps) && (count < maxiter))
  {
    new_params <- em_subfunction(FCST, OBS, PoC, B0, B1, C0, w, B_transform, tol.clip, count, CM2.iter)
    C0 <- new_params$C0
    w <- new_params$w
    error <- new_params$error
    count <- count + 1
  }

  lik <- get_log_lik(C0, w, OBS, FCST, B0, B1, PoC, B_transform, tol.clip)
  return(list(loglik=lik, w=w, C0=C0, count=count, max_error=max(abs(error))))
}

# ----------------------------------------------------------------------------------------
# EM helper functions

# either as equal weights or given starting values, ignoring forecast members which are unavailable
get_initial_weights <- function(start.w, avail) {
  if(is.null(start.w)){
    w <- sapply(avail, function(x) ifelse(x, 1/sum(avail), 0))
  }
  else{
    w <- sapply(seq_along(start.w), function(k) ifelse(avail[k], start.w[k]/sum(start.w[avail]), 0))
  }
  return(w)
}

em_subfunction <- function(FCST, OBS, PoC, B0, B1, C0, w, B_transform, tol.clip, count, CM2.iter) {
  ## E step
  z <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform, tol.clip)$z

  ## CM-1 step
  # new weights and variance deflation
  # n members = dim(FCST)[3]
  w.new <- sapply(1:dim(FCST)[3], function(k) {ifelse(any(!is.na(z[,,k])), mean(z[,,k], na.rm=TRUE), 0)})

  ## CM-1 step
  # Using the new w estimate, re-optimize C0
  if (count%%CM2.iter == 0) {
    opt <- optimize(get_log_lik, interval=c(0, 0.25), w=w.new, OBS=OBS, FCST=FCST, B0=B0, B1=B1, PoC=PoC, B_transform=B_transform, tol.clip=tol.clip,
                    maximum = T)
    C0.new <- opt$maximum
  } else C0.new <- C0

  # Complete list of changes to all w and c parameters
  error <- c(mapply("-", w.new, w), C0.new - C0)
  #-----------------------------

  return(list(C0=C0.new, error=error, w=w.new))
}

## E step subfunction
e_step <- function(w, C0, OBS, FCST, B0, B1, PoC, B_transform, tol.clip) {
  # Re-calculate beta density estimates based on current estimate for C0
  rhos <- array(mapply(get_rho, FCST, B0, B1, MoreArgs = list(B_transform=B_transform)), dim(FCST))
  gammas <- array(mapply(get_gamma, rhos, MoreArgs = list(C0=C0)), dim(FCST))
  # Using linear indexing, OBS gets recycled across members
  db <- array(mapply(dbeta_gamma_rho, OBS, gammas, rhos), dim(FCST))

  # z is an array with the first entry being day, second entry site, third entry forecast
  # ntime = dim(FCST)[1], nsite = dim(FCST)[2]
  z_num <-  array(mapply(get_z, OBS, PoC, db, rep(w, each=dim(FCST)[1]*dim(FCST)[2]), MoreArgs = list(tol.clip=tol.clip)), dim(FCST)) # Linear indexing, OBS is recycled across members, w explicitly expanded

  # sumz is weighted sum of density functions across the members. Rows are single training days, columns are sites
  sumz <- apply(z_num, MARGIN=c(1,2), FUN=function(z) ifelse(all(is.na(z)), NA, sum(z, na.rm=T))) # If all members are missing, return NA; else sum the others.
  z <- array(mapply("/", z_num, sumz), dim(z_num)) # sumz is recycled across members
  return(list(sumz=sumz, z=z))
}

# Define log likelihood subfunction for maximization step
get_log_lik <- function(C0, w, OBS, FCST, B0, B1, PoC, B_transform, tol.clip) {
  sumz <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform, tol.clip)$sumz
  return(sum(log(sumz), na.rm=T))
}

# Get z for single instance
# Returns NA for missing values and observations exactly at 0
# density is (PoC)*1[obs==1] + (1-PoC)*Beta(obs,a,b)*1[obs < 1]
# Uses tolerance to determine if obs == 1
get_z <- function(OBS, PoC, db, w, tol.clip) {
  if (is.na(OBS) | OBS==0) {return(NA)}
  else return(ifelse(abs(OBS-1) < tol.clip, w*PoC, w*(1-PoC)*db))
}

# Get PoC (probability of clipping) for single instance
# PoC is estimated from logit(PoC)=a0 + a1*f + a2*1[obs==1]
# This seems inverted from the PoP defintion from Will Kleiber's code, which reads to me like it should be 1-PoP
#' @param OBS single observation
#' @param FCST single forecast
#' @param A0 intercept of logistic regression of PoC
#' @param A1 slope of logistic regression of PoC
#' @param A2 coefficient of clipped forecast indicator [FCST==1] in logistic regression of PoC
get_poc <- function(FCST, A0, A1, A2, A_transform=NA) {
  if (is.na(FCST)) {return(NA)}
  else return(1/(1+exp(-(A0+A1*ifelse(typeof(A_transform)=='closure', A_transform(FCST), FCST)+A2*(FCST==1)))))
}

# Get rho parameter of beta distribution for single instance
#' @param FCST single forecast
#' @param B0 intercept of linear model of mean
#' @param B1 slope of linear model of mean
#' @param B_tranform Function of forecast transformation for B coefficients (optional), e.g. function(x) return(log(x)+1)
get_rho <- function(FCST, B0, B1, B_transform=NA) {
  if (is.na(FCST) | is.na(B0) | is.na(B1)) {return(NA)}
  # Estimate rho= mu (mean)
  mu <- B0 + B1*ifelse(typeof(B_transform)=='closure', B_transform(FCST), FCST)
  # Truncate mu just below 1, useful for those ensemble members that end up with positive B1
  if (mu >= 1)
      mu <- 1-1e-6
  return(mu)
}

# Get gamma parameter of beta distribution for single instance
#' @param mu beta mean value for this forecast (AKA rho)
#' @param C0 height parameter of quadratic model of variance
get_gamma <- function(mu, C0) {
  if (is.na(mu)) {return(NA)}

  # Estimate sigma (std deviation)
  sigma <- sqrt(-(C0/0.25)*(mu-0.5)^2 + C0)

  # Truncate at maximum theoretical value if necessary
  if (sigma >= sqrt(mu*(1-mu)))
    sigma <- sqrt(mu*(1-mu))-1e-6 #smidge less than the limit

  gamma <- mu*(1-mu)/sigma^2 - 1

  return(gamma)
}


# Get density of beta distribution with gamma, rho parameterization at given value
#' @param x Value or vector of values [0,1]
#' @param g Gamma = alpha + beta
#' @param rho rho = alpha/(alpha + beta)
dbeta_gamma_rho <- function(x, g, rho) {
  if (is.na(x) | is.na(g) | is.na(rho)) return(NA)
  gamma(g)/(gamma(g*rho)*gamma(g*(1-rho)))*x^(g*rho-1)*(1-x)^((1-rho)*g-1)
}

