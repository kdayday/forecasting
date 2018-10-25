# Methods for probabilistic forecast superclass
#------------------------------------------------------------------------------

#' Calculate number of sites in this spatial aggregate
#'
#' @param x A prob_forecast object
#' @return Number of sites
length.prob_forecast <- function(x){
  return(x$d)
}

#' Calculate forecast quantiles
#'
#' @param samples A column matrix of summed samples
#' @param quantile_density Numeric in (0,1), i.e., 0.1 to calculate quantiles at every 10%
#' @return A named numeric vector of estimated quantiles
calc_quantiles <- function(samples,
                     quantile_density=0.1) {
  if (quantile_density <= 0 | quantile_density >= 1) stop('Quantile density must be in (0,1).')
  quantiles <- stats::quantile(samples, probs=seq(0, 1 , quantile_density), type=1, names=TRUE)
  return(quantiles)
}

#' Calculate VaR and CVaR from sampled data. CVaR calculation is done directly from the samples, rather than estimated from a fitted distribution.
#'
#' @param samples Numeric vector
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return list of var, cvar
calc_cvar <- function(samples, epsilon=c(0.05, 0.95)) {
  if (any(epsilon <= 0) | any(epsilon >= 1)) stop("Epsilon's must be in (0,1).")
  var_low <- stats::quantile(samples, probs=epsilon[1], type=1, names=FALSE)
  cvar_low <- mean(samples[samples <= var_low])
  var_high <- stats::quantile(samples, probs=epsilon[2], type=1, names=FALSE)
  cvar_high <- mean(samples[samples >= var_high])
  return(list('cvar' = list('low'=cvar_low,'high' = cvar_high), 'var'=list('low'=var_low, 'high'=var_high)))
}

#' Calculate interval score for an interval from alpha/2 to 1-alpha/2. Negatively oriented
#' First term: sharpness, second and third terms: penalty for reliability.
#'
#' @param x A ts_forecast object
#' @param actual The realized value
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
calc_is <- function(x, actual, alpha) {
  if (alpha<=0 | alpha>=1) stop(paste('Alpha should be (0,1), given ', alpha, '.', sep=''))
  l <- unname(x$quantiles[paste(alpha/2*100, '%', sep='')])
  u <- unname(x$quantiles[paste((1-alpha/2)*100, '%', sep='')])
  if (is.na(l) | is.na(u)) stop("Requested quantile is not in the forecast's list of quantiles.")
  is <- (u-l) + (2/alpha)*(l-actual)*(actual < l) + (2/alpha)*(actual-u)*(actual > u)
}

#' Plot probabilistic forecast's estimated pdf (with kde)
#' Note that CVaR and VaR, while represented on the graph, are calculated directly from sampled data rather than estimated
#' from the kde results.
plot.prob_forecast <- function(x, ...) {
  # Assume data is power or irradiance and must be non-negative
  epdf <- stats::density(get_samples(x), from=0)
  plot(epdf, xlab='Power [W]', ylab='Probability',
       main='Estimated probability distribution', sub = paste("Location: ", x$location, ", Time:", x$time))

  # Color in tails above/below desired epsilon's
  i1 <- min(which(epdf$x >= x$var$low))
  i2 <- max(which(epdf$x <= x$var$high))
  graphics::lines(rep(x$var$low,times=2), c(0, epdf$y[i1]), col='black')
  graphics::polygon(c(0,epdf$x[1:i1],epdf$x[i1]), c(0, epdf$y[1:i1],0), col='red')
  graphics::text(epdf$x[i1], epdf$y[i1], paste("VaR: ", round(x$var$low,2), "\nCVaR: ", round(x$cvar$low,2)), pos=3)

  graphics::lines(rep(x$var$high,times=2), c(0, epdf$y[i2]), col='black')
  last <- length(epdf$x)
  graphics::polygon(c(epdf$x[i2], epdf$x[i2:last], epdf$x[last]), c(0, epdf$y[i2:last], 0), col='red')
  graphics::text(epdf$x[i2], epdf$y[i2], paste("VaR: ", round(x$var$high,2), "\nCVaR: ", round(x$cvar$high,2)), pos=3)
}

#' Check probabilistic forecast class
is.prob_forecast <- function(x) inherits(x, "prob_forecast")

#' Register generic sample function
#' @param x A prob_forecast object
get_samples <- function(x) {
    UseMethod("get_samples",x)
}

# Methods for aggregate probabilistic forecast class using vine copulas
#------------------------------------------------------------------------------


#' Initialize a probabilistic power forecast for a specific time point, using an n-dimensional vine copula.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#'
#' @param dat A matrix of training data [ntrain x nsites]
#' @param location A string
#' @param time A lubridate time stamp
#' @param training_transform_type One of "kde", "empirical" for transform of training data into uniform domain (default empirical)
#' @param results_transform_type One of "kde", "empirical" for transform of copula results back into variable domain (default kde)
#' @param n An integer, number of copula samples to take
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_vine_forecast <- function(dat, location, time,
                                  training_transform_type="empirical", results_transform_type='kde', n=3000, epsilon=c(0.05, 0.95)) {
  if (!is.numeric(n)) stop('n (number of samples) must be an integer.')
  if (dim(dat)[2] < 2) stop('Training data from more than 1 site required for vine copula forecast.')

  training_transforms <- apply(dat, MARGIN = 2, FUN=get_transform, transform_method=training_transform_type)
  # Results transforms must be subsequently updated if desired
  results_transforms <- apply(dat, MARGIN = 2, FUN=get_transform, transform_method=results_transform_type)

  uniform_dat <- mapply(function(n, t) {to_uniform(t, dat[,n])}, colnames(dat, do.NULL=FALSE), training_transforms)
  model <- rvinecopulib::vinecop(uniform_dat, family_set="all")

  # Initialize probabilistic forecast
  dat <- list(training_transforms = training_transforms,
              results_transforms = results_transforms,
              location = location,
              time = time,
              model = model,
              d = dim(dat)[2],
              n=n,
              epsilon=epsilon
              )
  x <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))

  # Complete probabilistic forecast by sampling and aggregating
  samples <- get_samples(x)
  x$quantiles <- calc_quantiles(samples)
  results <- calc_cvar(samples, epsilon)
  x$cvar <- results$cvar
  x$var<- results$var

  return(x)
}

#' Sample the vine copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @return A column matrix of aggregate powers
get_samples.prob_nd_vine_forecast <- function(x) {
  samples.u <- rvinecopulib::rvinecop(x$n, x$model)
  samples.xs <- matrix(nrow = x$n, ncol = length(x))
  for (i in 1:length(x)){
    samples.xs[,i] <- from_uniform(x$results_transforms, samples.u[,i])
  }
  samples.x <- rowSums(samples.xs)
  return(samples.x)
}

# Methods for aggregate probabilistic forecast class using Gaussian copulas
#------------------------------------------------------------------------------


#' Initialize a probabilistic power forecast for a specific time point, using an n-dimensional Gaussian copula.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#'
#' @param dat A matrix of training data, [ntrain x 1] for sites, [ntrain x nsites] for regional or total forecasts
#' @param location A string
#' @param time A lubridate time stamp
#' @param n An integer, number of copula samples to take
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_gaussian_forecast <- function(dat, location, time,
                                      n=3000, epsilon=c(0.05, 0.95)) {
  stop('Not implemented')
}

#' Sample the gaussian copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @param x A prob_forecast object
#' @return A column matrix of aggregate powers
get_samples.prob_nd_gaussian_forecast <- function(x) {
  stop('Not implemented')
}

# Methods for aggregate probabilistic forecast class using empirical copulas
#------------------------------------------------------------------------------


#' Initialize a probabilistic power forecast for a specific time point, using an n-dimensional empirical copula.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#'
#' @param dat A matrix of training data, [ntrain x nsites]
#' @param location A string
#' @param time A lubridate time stamp
#' @param n An integer, number of copula samples to take
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_empirical_forecast <- function(dat, location, time,
                                      n=3000, epsilon=c(0.05, 0.95)) {
  if (dim(dat)[2] < 2) stop('Training data from more than 1 site required for empirical copula forecast.')
  stop('Not implemented')
}

#' Sample the empirical copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @param x A prob_forecast object
#' @return A column matrix of aggregate powers
get_samples.prob_nd_empirical_forecast <- function(x) {
  stop('Not implemented')
}

# Methods for univariate probabilistic forecast class
#------------------------------------------------------------------------------


#' Initialize a univariate probabilistic power forecast for a specific time point.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#'
#' @param dat A matrix of training data, [ntrain x 1]
#' @param location A string
#' @param time A lubridate time stamp
#' @param n An integer, number of copula samples to take
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_1d_site_forecast <- function(dat, location, time,
                                       n=3000, epsilon=c(0.05, 0.95)) {
  if (dim(dat)[2] > 1) stop('Training data must be of dimensions [ntrain x 1] for univariate forecasts.')
  stop('Not implemented')
}

#' Sample the empirical copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @param x A prob_forecast object
#' @return A column matrix of aggregate powers
get_samples.prob_1d_site_forecast <- function(x) {
  stop('Not implemented')
}
