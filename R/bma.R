

# Fit BMA coefficients for beta distributions with discrete component at 1 for each ensemble member
# Data should be pre-processed and normalized to [0,1] so that clipped values are exactly 1
#' @param tel Vector of training telemetry data on [0,1]
#' @param ens Matrix of training ensemble member data [time x member] on [0,1]
#' @param lr_formula Formula in terms of x,y for logistic regression model, defaults to "y ~ x". Requires a negative x intercept to model PoC < 0.5.
#' @param A_transform A function for transforming forecast data before logistic regression to get a's (optional)
#' @param lm_formula Formula in terms of x,y for linear regression model, defaults to "y ~ x + 0"
#' @param B_transform A function for transforming forecast data before linear regression to get b's (optional)
#' @param percent_clipping_threshold [0,1] Power is designated as clipped when above this percentage of the maximum power
#' @param tol A tolerance for determining if normalized values fall are <=1 (defaults to 0.001)
#' @param ... Optional arguments to pass to EM algorithm
#' @return A list for a discrete-continuous mixture model with beta distribution
beta1_ens_models <- function(tel, ens, lr_formula= y ~ x, A_transform=NA, lm_formula= y ~ x + 0, B_transform=NA, percent_clipping_threshold=0.995, tol=1e-3, ...) {
  if (any(tel < 0, na.rm=T) | any(tel - 1>tol, na.rm=T)) stop('Telemetry must be normalized to [0,1] to apply beta model.')
  if (any(ens < 0, na.rm=T) | any(ens - 1>tol , na.rm=T)) stop('All forecasts must be normalized to [0,1] to apply beta model.')
  if (percent_clipping_threshold < 0.95  | percent_clipping_threshold > 1) stop('Percent clipping threshold must be <=1. Much greater than 0.95 is recommended.')
  if (length(tel) != dim(ens)[1]) stop("Must have same number of telemetry and forecast time-points.")

  # 1. Logistic regression for a's
  # Returns a list of model summaries, one per member
  mem_discrete_models <-lapply(seq_len(dim(ens)[2]), function(i) get_lr(ens[,i], tel=tel, form=lr_formula, A_transform = A_transform,
                                                                        percent_clipping_threshold=percent_clipping_threshold))
  A0 <- sapply(mem_discrete_models, function(m) return(unname(m$coefficients["(Intercept)"])))
  A1 <- sapply(mem_discrete_models, function(m) return(unname(m$coefficients["x"])))

  # 2. linear regression for b's.
  mem_mean_models <- lapply(seq_len(dim(ens)[2]), function(i) get_lm(ens[,i], tel=tel, form=lm_formula, B_transform = B_transform,
                                                                     percent_clipping_threshold=percent_clipping_threshold))
  # Force intercept to 0 if it is not included in the model
  B0 <- sapply(mem_mean_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Estimate"], 0)))
  # Force slope NA if it is missing -- this can happen when the member training data is always 0, returning in singularities in the glm fit.
  B1 <- sapply(mem_mean_models, function(m) return(ifelse("x" %in% rownames(m$coefficients), m$coefficients["x", "Estimate"], NA)))

  fit_statistics <- data.frame("A0 p-value"=sapply(mem_discrete_models, function(m) return(unname(m$prob["(Intercept)"]))),
                               "A1 p-value"=sapply(mem_discrete_models, function(m) return(unname(m$prob["x"]))),
                               "A AIC"=sapply(mem_discrete_models, function(m) return(m$aic)),
                               "B0 p-value"=sapply(mem_mean_models, function(m) return(ifelse("(Intercept)" %in% rownames(m$coefficients), m$coefficients["(Intercept)", "Pr(>|t|)"], NA))),
                               "B1 p-value"=sapply(mem_mean_models, function(m) return(ifelse("x" %in% rownames(m$coefficients), m$coefficients["x", "Pr(>|t|)"], NA))),
                               "B r-squared"=sapply(mem_mean_models, function(m) return(m$r.squared))) # If singularities occur in fitting, r.squared value with be 0.

  # 3. ME algorithm for w's and c's
  # In future, the above could be expanded to have unique site values. For now, matrices are expanded to arrays assuming a "single" site (i.e., global values).
  ntime <- dim(ens)[1]
  nens <- dim(ens)[2]
  array_dims <- c(ntime, 1, nens)
  tmp <- em(FCST=array(ens, dim=array_dims), OBS=array(tel, dim=c(ntime, 1)), A0=array(rep(A0, each=ntime), dim=array_dims),
            A1=array(rep(A1, each=ntime), dim=array_dims), A2=0, B0=array(rep(B0, each=ntime), dim=array_dims),
            B1=array(rep(B1, each=ntime), dim=array_dims), A_transform=A_transform, B_transform=B_transform, percent_clipping_threshold=percent_clipping_threshold, ...)

  return(list(A0=A0, A1=A1, A2=0, B0=B0, B1=B1, C0=tmp$C0, w=tmp$w, fit_statistics=fit_statistics, log_lik=tmp$log_lik,
              A_transform=A_transform, B_transform=B_transform, em_count=tmp$count, em_error=tmp$max_error, percent_clipping_threshold=percent_clipping_threshold))
}


# Do logistic regression to get a single ensemble member's 'a' coefficients.
#' Uses a penalized/Firth linear regression to handle quasi- or complete separation
# Telemetry within given tolerance of 1 are assumed to be clipped
# glm and lm omit NA values by default
#' @param fc Vector of training forecast data on [0,1]
#' @param tel Vector of telemetry data on [0,1]
#' @param form Formula in terms of x,y for logistic regression model
#' @param A_transform A function for pre-transforming the forecast data (optional), e.g. function(X) return(log(x)+1)
#' @param percent_clipping_threshold Percentage threshold for determining if clipping is occurring
#' @return Summary of the glm model
get_lr <- function(fc, tel, form, A_transform, percent_clipping_threshold){
  clipped <- as.integer(tel >= percent_clipping_threshold)
  if (typeof(A_transform)=="closure") fc <- sapply(fc, FUN = A_transform)
  # Check for complete separation (i.e., there are no clipped data or all clipped data in the training set)
  if (all(clipped[!is.na(tel) & !is.na(fc) & fc != 0]==0)) {
    # If no data points in the training set are clipped, set logistic regression intercept to -Inf
    coefficients <- c(-Inf, 0)
    names(coefficients) <- c("(Intercept)", "x")
    probs <- c(NA, NA)
    names(probs) <- c("(Intercept)", "x")
    return(list(aic=NA, coefficients=coefficients, prob=probs))
  } else if (all(clipped[!is.na(tel) & !is.na(fc) & fc != 0]==1)) {
    # If all data points in the training set are clipped, set logistic regression intercept to +Inf
    coefficients <- c(Inf, 0)
    names(coefficients) <- c("(Intercept)", "x")
    probs <- c(NA, NA)
    names(probs) <- c("(Intercept)", "x")
    return(list(aic=NA, coefficients=coefficients, prob=probs))
  } else {
    fit <- logistf::logistf(form, data=data.frame(x=fc, y=clipped))
    return(c(fit, aic=extractAIC(fit)[2]))
  }
}


# Get the linear model for a single ensemble member's 'b' coefficients
# Telemetry within given tolerance of 1 are assumed to be clipped
# glm and lm omit NA values by default
#' @param fc Vector of training forecast data on [0,1]
#' @param tel Vector of telemetry data on [0,1]
#' @param form Formula in terms of x,y for logistic regression model
#' @param B_transform A function for pre-transforming the forecast data (optional), e.g. function(X) return(log(x)+1)
#' @param percent_clipping_threshold Percentage threshold for determining if clipping is occurring
#' @return Summary of the lm model
get_lm <- function(fc, tel, form, B_transform, percent_clipping_threshold){
  unclipped <- tel < percent_clipping_threshold
  if (typeof(B_transform)=="closure") fc <- sapply(fc, FUN = B_transform)
  return(summary(lm(form, data=data.frame(x=fc[unclipped], y=tel[unclipped]))))
}

# Expectation-maximization function, modified from code courtesy of Will Kleiber
em <- function(FCST, OBS, A0, A1, A2, B0, B1, A_transform, B_transform, percent_clipping_threshold, C0=0.06, eps=1e-005, maxiter=1000, CM2.iter=50, start.w=NULL)

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
    new_params <- em_subfunction(FCST, OBS, PoC, B0, B1, C0, w, B_transform, percent_clipping_threshold, count, CM2.iter)
    C0 <- new_params$C0
    w <- new_params$w
    error <- new_params$error
    count <- count + 1
  }

  lik <- get_log_lik(C0, w, OBS, FCST, B0, B1, PoC, B_transform, percent_clipping_threshold)
  return(list(log_lik=lik, w=w, C0=C0, count=count, max_error=max(abs(error))))
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

em_subfunction <- function(FCST, OBS, PoC, B0, B1, C0, w, B_transform, percent_clipping_threshold, count, CM2.iter) {
  ## E step
  z <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform, percent_clipping_threshold)$z

  ## CM-1 step
  # new weights and variance deflation
  # n members = dim(FCST)[3]
  w.new <- sapply(1:dim(FCST)[3], function(k) {ifelse(any(!is.na(z[,,k])), mean(z[,,k], na.rm=TRUE), 0)})

  ## CM-2 step
  # Using the new w estimate, re-optimize C0
  if (count%%CM2.iter == 0) {
    opt <- optimize(get_log_lik, interval=c(0, 0.25), w=w.new, OBS=OBS, FCST=FCST, B0=B0, B1=B1, PoC=PoC, B_transform=B_transform,
                    percent_clipping_threshold=percent_clipping_threshold, maximum = T)
    C0.new <- opt$maximum
  } else C0.new <- C0

  # Complete list of changes to all w and c parameters
  error <- c(mapply("-", w.new, w), C0.new - C0)
  #-----------------------------

  return(list(C0=C0.new, error=error, w=w.new))
}

## E step subfunction
e_step <- function(w, C0, OBS, FCST, B0, B1, PoC, B_transform, percent_clipping_threshold) {
  # Re-calculate beta density estimates based on current estimate for C0
  rhos <- array(mapply(get_rho, FCST, B0, B1, MoreArgs = list(B_transform=B_transform)), dim(FCST))
  gammas <- array(mapply(get_gamma, rhos, MoreArgs = list(C0=C0)), dim(FCST))

  # z is an array with the first entry being day, second entry site, third entry forecast
  # ntime = dim(FCST)[1], nsite = dim(FCST)[2]
  # Linear indexing, OBS is recycled across members, w explicitly expanded
  z_num <-  array(mapply(get_weighted_probability, OBS, PoC, gammas, rhos, rep(w, each=dim(FCST)[1]*dim(FCST)[2]), MoreArgs = list(percent_clipping_threshold=percent_clipping_threshold)), dim(FCST))

  # sumz is weighted sum of density functions across the members. Rows are single training days, columns are sites
  sumz <- apply(z_num, MARGIN=c(1,2), FUN=function(z) ifelse(all(is.na(z)), NA, sum(z, na.rm=T))) # If all members are missing, return NA; else sum the others.
  z <- array(mapply("/", z_num, sumz), dim(z_num)) # sumz is recycled across members
  return(list(sumz=sumz, z=z))
}

# Define log likelihood subfunction for maximization step
get_log_lik <- function(C0, w, OBS, FCST, B0, B1, PoC, B_transform, percent_clipping_threshold) {
  sumz <- e_step(w, C0, OBS, FCST, B0, B1, PoC, B_transform, percent_clipping_threshold)$sumz
  return(sum(log(sumz), na.rm=T))
}

# Get z=weighted probability for single instance
# Returns NA for missing values and observations exactly at 0
# Defines clipping with a percentage threshold, lambda
# density is (PoC/(1-lambda))*1[obs>=lambda] + ((1-PoC)/CDF(lambda))*Beta(obs,a,b)*1[obs < lambda]
get_weighted_probability <- function(OBS, PoC, gamma, rho, w, percent_clipping_threshold) {
  if (is.na(OBS) | OBS==0) {return(NA)}
  if (OBS >= percent_clipping_threshold) return(w*PoC/(1-percent_clipping_threshold))
  else {
    alpha <- rho * gamma
    beta <- gamma * (1-rho)
    db <- stats::dbeta(OBS, alpha, beta)
    threshold_CD <- stats::pbeta(percent_clipping_threshold, alpha, beta) # Scaling factor: cumulative density at the threshold
    return((w*(1-PoC)/threshold_CD)*db)
  }
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
  if (mu <= 0 )
      mu <- 1e-6
  return(mu)
}

# Get gamma parameter of beta distribution for single instance
#' @param mu beta mean value for this forecast (AKA rho)
#' @param C0 height parameter of quadratic model of variance
get_gamma <- function(mu, C0) {
  # Return NA if mu is missing or on the boundary
  if (is.na(mu) | mu==0) {return(NA)}

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
  if (all(is.na(x)) | is.na(g) | is.na(rho)) return(NA)
  # This is calculated using the log gamma function instead of the pure gamma, to avoid intermediate Inf values near the domain extremes.
  exp(lgamma(g)-(lgamma(g*rho)+lgamma(g*(1-rho))))*x^(g*rho-1)*(1-x)^((1-rho)*g-1)
}

