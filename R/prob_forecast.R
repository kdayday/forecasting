# Methods for probabilistic forecast superclass
#------------------------------------------------------------------------------

#' Calculate number of sites in this spatial aggregate
#'
#' @param x A prob_forecast object
#' @return Number of sites
length.prob_forecast <- function(x){
  return(x$d)
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



#' Plot probabilistic forecast's quantiles
plot.prob_forecast <- function(x) {
  plot(x$quantiles, seq(0, 1, length.out = length(x$quantiles)), xlab='Power [W]', ylab='Cumulative Density',
       main='Quantiles', sub = paste("Location: ", x$location, ", Time:", x$time))
}

#' Check probabilistic forecast class
is.prob_forecast <- function(x) inherits(x, "prob_forecast")

#' Register generic sample function
#' @param x A prob_forecast object
get_1d_samples <- function(x) {
    UseMethod("get_1d_samples",x)
}

#' Register generic joint density function
#' @param x A Gaussian or vine copula n-dimension prob_forecast object
get_joint_density_grid <- function(x, ...) {
  UseMethod("get_joint_density_grid",x)
}

#' Register generic quantiles function
#' @param x A prob_forecast object
calc_quantiles <- function(x, ...) {
  UseMethod("calc_quantiles",x)
}

#' Register generic VaR/CVaR function
#' @param x A prob_forecast object
calc_cvar <- function(x, ...) {
  UseMethod("calc_cvar",x)
}


#' Register generic plot function
#' @param x A prob_forecast object
plot_pdf <- function(x, ...) {
  UseMethod("plot_pdf",x)
}

# Methods for aggregate probabilistic forecast class using vine copulas
#------------------------------------------------------------------------------


#' Initialize a probabilistic power forecast for a specific time point, using an n-dimensional vine copula.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#'
#' @param dat A matrix of training data [ntrain x nsites]
#' @param location A string
#' @param time A lubridate time stamp
#' @param training_transform_type Transform of training data into uniform domain (see marg_transform "method")
#' @param results_transform_type Transform of copula results back into variable domain (see marg_transform "method")
#' @param n An integer, number of copula samples to take
#' @param ... optional arguments to the marginal estimator
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_vine_forecast <- function(dat, location, time,
                                  training_transform_type="empirical", results_transform_type='empirical', n=3000, ...) {
  if (!is.numeric(n)) stop('n (number of samples) must be an integer.')
  if (class(dat)!='matrix') stop('Input data must be a matrix')
  if (dim(dat)[2] < 2) stop('Training data from more than 1 site required for vine copula forecast.')

  tr <- calc_transforms(dat, training_transform_type, results_transform_type, ...)
  training_transforms <- tr$training
  results_transforms <- tr$results

  uniform_dat <- mapply(function(n, t) {to_uniform(t, dat[,n])}, colnames(dat, do.NULL=FALSE), training_transforms)
  model <- rvinecopulib::vinecop(uniform_dat, family_set="all")

  # Initialize probabilistic forecast
  dat <- list(training_transforms = training_transforms,
              results_transforms = results_transforms,
              location = location,
              time = time,
              model = model,
              d = dim(dat)[2],
              n=n
              )
  x <- structure(dat, class = c("prob_forecast", "prob_nd_vine_forecast"))

  # Complete probabilistic forecast by sampling and aggregating
  x$quantiles <- calc_quantiles(x)

  return(x)
}

#' Check class
is.prob_nd_vine_forecast <- function(x) inherits(x, "prob_nd_vine_forecast")

#' Calculate lists of variable-to-uniform domain transforms for all dimensions
#'
#' @param dat training data matrix
#' @param training_transform_type Transform of training data into uniform domain (see marg_transform "method")
#' @param results_transform_type Transform of copula results back into variable domain (see marg_transform "method")
#' @param ... Optional arguments to marg_transform
#' @return list of "training" and "results" transforms to use.
calc_transforms <- function(dat, training_transform_type, results_transform_type, ...) {
  training <- lapply(seq_len(dim(dat)[2]), FUN=get_transform_with_unique_xmin_max, dat=dat, method=training_transform_type, ...)
  # Results transforms must be subsequently updated if desired
  if (results_transform_type==training_transform_type) {
    results <- training
  } else {
    results <- lapply(seq_len(dim(dat)[2]), FUN=get_transform_with_unique_xmin_max, dat=dat, method=results_transform_type, ...)
  }
  return(list('training'=training, 'results'=results))
}

#' Subfunction for calc_transforms to cycle thorugh xmin and xmax if they are given uniquely for each dimension
#'
#' @param idx Column index of dat
#' @param dat training data matrix over all the dimensions
#' @param method marg_transform method
#' @param ... Optional arguments to marg_transform, including potentially xmin or xmax in either scalar or vector form
get_transform_with_unique_xmin_max <- function(idx, dat, method, ...) {
  args <- list(...)
  # Use unique xmin/xmax values if vectors are given
  if ('xmin' %in% names(args) & length(args[['xmin']]) > 1) {args[['xmin']] <- args[['xmin']][idx]}
  if ('xmax' %in% names(args) & length(args[['xmax']]) > 1) {args[['xmax']] <- args[['xmax']][idx]}
  return(do.call(marg_transform, c(list(dat[,idx], method), args))) # repackage arguments into single list for do.call
}

#' Sample the vine copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @return A column matrix of aggregate powers
get_1d_samples.prob_nd_vine_forecast <- function(x) {
  samples.u <- rvinecopulib::rvinecop(x$n, x$model)
  samples.xs <- matrix(nrow = x$n, ncol = length(x))
  for (i in 1:length(x)){
    samples.xs[,i] <- from_uniform(x$results_transforms[[i]], samples.u[,i])
  }
  samples.x <- rowSums(samples.xs)
  return(samples.x)
}

#' Calculate forecast quantiles from samples of the vine copula
#'
#' @param x prob_nd_vine_forecast object
#' @param samples (optional) previously obtained samples to use instead of new sampling, e.g. for coordination with cVaR calculation
#' @param quantile_density Numeric in (0,1), i.e., 0.1 to calculate quantiles at every 10%
#' @return A named numeric vector of estimated quantiles
calc_quantiles.prob_nd_vine_forecast <- function(x, samples=NA, quantile_density=0.1) {
  if (quantile_density <= 0 | quantile_density >= 1) stop('Bad input. Quantile density must be in (0,1).')

  if (!(is.numeric(samples))) {samples <- get_1d_samples(x)}
  quantiles <- stats::quantile(samples, probs=seq(0, 1 , quantile_density), type=1, names=TRUE)
  return(quantiles)
}

#' Calculate VaR and CVaR from sampled data. CVaR calculation is done directly from the samples, rather than estimated from a fitted distribution.
#'
#' @param x prob_nd_vine_forecast object
#' @param samples (optional) previously obtained samples to use instead of new sampling, e.g. for coordination with quantiles calculation
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return list of var, cvar
calc_cvar.prob_nd_vine_forecast <- function(x, samples=NA, epsilon=c(0.05, 0.95)) {
  if (!(is.numeric(samples))) {samples <- get_1d_samples(x)}
  if (any(epsilon <= 0) | any(epsilon >= 1)) stop("Bad input. Epsilon's must be in (0,1).")

  var_low <- stats::quantile(samples, probs=epsilon[1], type=1, names=FALSE)
  cvar_low <- mean(samples[samples <= var_low])
  var_high <- stats::quantile(samples, probs=epsilon[2], type=1, names=FALSE)
  cvar_high <- mean(samples[samples >= var_high])
  return(list('cvar' = list('low'=cvar_low,'high' = cvar_high), 'var'=list('low'=var_low, 'high'=var_high)))
}

#' Estimate the joint probability density
#'
#' @param x A prob_nd_vine_forecast object
#' @param k Integer or vector of number of samples to take along each dimension. If given an integer, all dimensions have the same number of samples.
#' @return an joint probability density array
get_joint_density_grid.prob_nd_vine_forecast <- function(x, k=100) {
  # Get evaluation grid
  eval_points <- get_variable_domain_grid(x, k)

  # Calculate copula density
  eval_points_copula <- mapply(FUN=function(c, n) {to_uniform(c, eval_points[,n])}, x$results_transforms,
                               colnames(eval_points, do.NULL=FALSE))
  copula_density <- rvinecopulib::dvinecop(eval_points_copula, x$model)

  # Calculate the marginal densities
  dmarg <- mapply(FUN=function(c, n) {to_probability(c, eval_points[,n])}, x$results_transforms,
                  colnames(eval_points, do.NULL=FALSE))

  # Calculate joint density as product the copula density and marginal densities
  pdf <- copula_density * apply(dmarg, 1 , prod) # multiply across rows
  return(list('grid_points'=eval_points, 'd'=pdf))
}

#' Get a grid of evaluation points in the variable domain, based on the range of the marginal distribution estimation
#'
#' @param x A prob_nd_vine_forecast object
#' @param k Integer or vector of number of samples to take along each dimension. If given an integer, all dimensions have the same number of samples.
#' @return a grid of points size k^n
get_variable_domain_grid <- function(x, k) {
  if (length(k)==1){
    k <- rep(k, times=length(x))
  } else if (length(k) != length(x)) stop('Bad input. k must be length 1 or of same dimension as the forecast')
  if (any(k < 2)) stop('Bad input. All values in k must be at least 2.')

  grid_list <-mapply(x$results_transforms, k, FUN=function(tr, k) {seq(tr$xmin, tr$xmax, length.out=k)}, SIMPLIFY = FALSE)
  pts <- expand.grid(grid_list)
  return(pts)
}

#' Plot probabilistic forecast's estimated pdf (with kde)
#' Note that CVaR and VaR, while represented on the graph, are calculated directly from sampled data rather than estimated
#' from the kde results.
#'
#' @param x prob_nd_vine_forecast object
plot_pdf.prob_nd_vine_forecast <- function(x, cvar=FALSE, epsilon=c(0.05, 0.95)) {
  # Assume data is power or irradiance and must be non-negative
  samples <- get_1d_samples(x)
  epdf <- stats::density(samples, from=0)
  plot(epdf, xlab='Power [W]', ylab='Probability',
       main='Estimated probability density', sub = paste("Location: ", x$location, ", Time:", x$time))

  if (cvar){# Color in tails above/below desired epsilon's
    cvar_info <- calc_cvar(x, samples=samples, epsilon=epsilon)
    i1 <- min(which(epdf$x >= cvar_info$var$low))
    i2 <- max(which(epdf$x <= cvar_info$var$high))
    graphics::lines(rep(cvar_info$var$low,times=2), c(0, epdf$y[i1]), col='black')
    graphics::polygon(c(0,epdf$x[1:i1],epdf$x[i1]), c(0, epdf$y[1:i1],0), col='red')
    graphics::text(epdf$x[i1], epdf$y[i1], paste("VaR: ", round(cvar_info$var$low,2), "\nCVaR: ", round(cvar_info$cvar$low,2)), pos=3)

    graphics::lines(rep(cvar_info$var$high,times=2), c(0, epdf$y[i2]), col='black')
    last <- length(epdf$x)
    graphics::polygon(c(epdf$x[i2], epdf$x[i2:last], epdf$x[last]), c(0, epdf$y[i2:last], 0), col='red')
    graphics::text(epdf$x[i2], epdf$y[i2], paste("VaR: ", round(cvar_info$var$high,2), "\nCVaR: ", round(cvar_info$cvar$high,2)), pos=3)
  }
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
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_gaussian_forecast <- function(dat, location, time, n=3000, ...) {
  stop('Not implemented')
}

#' Check class
is.prob_nd_gaussian_forecast <- function(x) inherits(x, "prob_nd_gaussian_forecast")

#' Sample the gaussian copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @param x A prob_forecast object
#' @return A column matrix of aggregate powers
get_1d_samples.prob_nd_gaussian_forecast <- function(x) {
  stop('Not implemented')
}


#' Estimate the joint probability density
#'
#' @param x A prob_nd_gaussian_forecast object
#' @param k Integer or vector of number of samples to take along each dimension. If given an integer, all dimensions have the same number of samples.
#' @return an joint probability density array
get_joint_density_grid.prob_nd_gaussian_forecast <- function(x, k=100) {
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
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_nd_empirical_forecast <- function(dat, location, time, n=3000) {
  if (dim(dat)[2] < 2) stop('Training data from more than 1 site required for empirical copula forecast.')
  stop('Not implemented')
}

#' Check class
is.prob_nd_empirical_forecast <- function(x) inherits(x, "prob_nd_empirical_forecast")


#' Sample the empirical copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @param x A prob_forecast object
#' @return A column matrix of aggregate powers
get_1d_samples.prob_nd_empirical_forecast <- function(x) {
  stop('Not implemented')
}

# Methods for univariate probabilistic forecast class
#------------------------------------------------------------------------------

#' Initialize a univariate probabilistic power forecast for a specific time point by ranking ensemble members.
#' Its rank quantiles are the basic quantiles estimated from the ensemble members; the quantiles are the rank quantiles iterpolated
#' to the quantiles of interest (i.e., 10%, 20%, etc...)
#'
#' @param dat A numeric vector of ensemble members
#' @param location A string
#' @param time A lubridate time stamp
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_1d_rank_forecast <- function(dat, location, time, ...) {
  if (!is.vector(dat)) stop("Input data must be a vector.")

  # Initialize probabilistic forecast
  dat <- list(location = location,
              rank_quantiles = list(x=sort(dat), y=(rank(sort(dat))-1)/(length(dat)-1)),
              time = time,
              d = 1)

  x <- structure(dat, class = c("prob_forecast", "prob_1d_rank_forecast"))
  x$quantiles <- calc_quantiles(x)

  return(x)
}

#' Check class
is.prob_1d_rank_forecast <- function(x) inherits(x, "prob_1d_rank_forecast")

#' Calculate forecast quantiles
#'
#' @param x prob_1d_rank_forecast object
#' @param quantile_density Numeric in (0,1), i.e., 0.1 to calculate quantiles at every 10%
#' @return A named numeric vector of estimated quantiles
calc_quantiles.prob_1d_rank_forecast <- function(x, quantile_density=0.1) {
  if (quantile_density <= 0 | quantile_density >= 1) stop('Bad input. Quantile density must be in (0,1).')
  yseq <- seq(0, 1, by=quantile_density)
  xseq <- stats::approx(x=x$rank_quantiles$y,  y=x$rank_quantiles$x, xout=yseq)$y

  names(xseq) <- sapply(yseq, FUN=function(y) return(paste(y*100, "%", sep='')))
  return(xseq)
}


#' Not implemented
#'
#' @param x prob_1d_rank_forecast object
plot_pdf.prob_1d_rank_forecast <- function(x) {
  stop("Not implemented for 1D rank forecasts.")
}

# ---------------------------------------------------------------------------------------------

#' Initialize a univariate probabilistic power forecast for a specific time point using kernel density estimation. See kde_methods.R for more details
#'
#' @param dat A vector of ensemble members
#' @param location A string
#' @param time A lubridate time stamp
#' @param method KDE method selection, see kde_methods.R for details
#' @param ... Additional parameters passed on the KDE method
#' @return An n-dimensional probabilistic forecast object from vine copulas
prob_1d_kde_forecast <- function(dat, location, time, method='geenens', ...) {
  if (!is.vector(dat)) stop("Input data must be a vector.")

  func <- kde_lookup(method)
  # Get selected KDE
  model <- func(dat, ...)

  # Initialize probabilistic forecast
  dat <- list(location = location,
              time = time,
              d = 1,
              model=model
              )
  x <- structure(dat, class = c("prob_forecast", "prob_1d_kde_forecast"))

  # Complete probabilistic forecast by sampling and aggregating
  x$quantiles <- calc_quantiles(x)
  return(x)
}

#' Check class
is.prob_1d_kde_forecast <- function(x) inherits(x, "prob_1d_kde_forecast")

#' Calculate forecast quantiles from KDE
#'
#' @param x prob_1d_kde_forecast object
#' @param quantile_density Numeric in (0,1), i.e., 0.1 to calculate quantiles at every 10%
#' @return A named numeric vector of estimated quantiles
calc_quantiles.prob_1d_kde_forecast <- function(x, quantile_density=0.1) {
  if (quantile_density <= 0 | quantile_density >= 1) stop('Bad input. Quantile density must be in (0,1).')
  yseq <- seq(0, 1, by=quantile_density)
  xseq <- stats::approx(x=x$model$u,  y=x$model$x, xout=yseq)$y

  names(xseq) <- sapply(yseq, FUN=function(y) return(paste(y*100, "%", sep='')))
  return(xseq)
}

#' Calculate VaR and CVaR by trapezoidal integration.
#'
#' @param x prob_1d_kde_forecast object
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return list of var, cvar
calc_cvar.prob_1d_kde_forecast <- function(x, epsilon=c(0.05, 0.95)) {
  if (any(epsilon <= 0) | any(epsilon >= 1)) stop("Bad input. Epsilon's must be in (0,1).")

  var_low <- stats::approx(x=x$model$u, y=x$model$x, xout=epsilon[1])$y
  idx <- min(which(x$model$x >= var_low))
  cvar_low <- (1/epsilon[1])*pracma::trapz(x$model$x[1:idx]*x$model$d[1:idx], x$model$x[1:idx])

  var_high <- stats::approx(x=x$model$u, y=x$model$x, xout=epsilon[2])$y
  idx <- max(which(x$model$x <= var_high))
  cvar_high <- (1/(1-epsilon[2]))*pracma::trapz(x$model$x[-(1:(idx-1))]*x$model$d[-(1:(idx-1))], x$model$x[-(1:(idx-1))])

  return(list('cvar' = list('low'=cvar_low,'high' = cvar_high), 'var'=list('low'=var_low, 'high'=var_high)))
}


#' Plot probability density
#'
#' @param x prob_1d_kde_forecast object
plot_pdf.prob_1d_kde_forecast <- function(x) {
  plot(fc_kde$model$x, fc_kde$model$d, xlab='Power [W]', ylab='Probability density', sub = paste("Location: ", x$location, ", Time:", x$time),
       type='l', lwd=2)
}
