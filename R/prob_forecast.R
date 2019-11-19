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
#' @param x A prob_forecast object
#' @param actual The realized value
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
IS <- function(x, tel, alpha) {
  eps <- 1e-6
  if (alpha<=0 | alpha>=1) stop(paste('Alpha should be (0,1), given ', alpha, '.', sep=''))
  l <- x$quantiles$x[abs(x$quantiles$q-alpha/2)<eps]
  u <- x$quantiles$x[abs(x$quantiles$q-(1-alpha/2))<eps]
  if (length(l)==0 | length(u)==0) stop("Requested quantile is not in the forecast's list of quantiles.")
  is <- (u-l) + (2/alpha)*(l-tel)*(tel < l) + (2/alpha)*(tel-u)*(tel > u)
}

#' Calculate sharpness for an interval from alpha/2 to 1-alpha/2. Negatively oriented
#'
#' @param x A prob_forecast object
#' @param tel (unneeded) For calling signature consistency with other metrics
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
sharpness <- function(x, tel, alpha) {
  eps <- 1e-6
  if (alpha<=0 | alpha>=1) stop(paste('Alpha should be (0,1), given ', alpha, '.', sep=''))
  l <- x$quantiles$x[abs(x$quantiles$q-alpha/2)<eps]
  u <- x$quantiles$x[abs(x$quantiles$q-(1-alpha/2))<eps]
  if (length(l)==0 | length(u)==0) stop("Requested quantile is not in the forecast's list of quantiles.")
  sharpness <- u-l
}

#' Estimate CRPS
#'
#' @param x A prob_forecast object
#' @param tel Telemetry value
CRPS <- function(x, tel) {
  y <- x$quantiles$x

  # CRPS broken down into two parts below and above tel to simplify Heaviside evaluation,
  # plus one of two optional rectangles representing additional CRPS area if the telemetry was an outlier outside the forecasted quantiles
  crps <- ifelse(tel < min(y), min(y)-tel, 0) + # lower outlier rectangular area of height 1 (implied)
          pracma::trapz(y[y <= tel], x$quantiles$q[y <= tel]^2) + # Area below heaviside step
          pracma::trapz(y[y >= tel], (x$quantiles$q[y >= tel]-1)^2) + # Area above heaviside step
          ifelse(tel > max(y), tel-max(y), 0)# upper outlier rectangular area of height 1 (implied)
}

#' Plot probabilistic forecast's quantiles
plot.prob_forecast <- function(x) {
  plot(x$quantiles$x, x$quantiles$q, xlab='Power [MW]', ylab='Cumulative Distribution',
      sub = paste("Location: ", x$location, ", Time:", x$time))
}

#' Plot probabilistic forecast's CRPS illustration
plot_crps <- function(x, tel) {
  if (!is.prob_forecast(x)) stop("x must be a prob_forecast object to plot CRPS")
  y <- x$quantiles$x

  if (tel < min(y)) {
    upper_area.x <- c(tel, min(y), x$quantiles$x[y >= tel])
    upper_area.y <- c(1, 1, (x$quantiles$q[y >= tel]-1)^2)
  } else {
    upper_area.x <- x$quantiles$x[y >= tel]
    upper_area.y <- (x$quantiles$q[y >= tel]-1)^2
  }
  if (tel > max(y)) {
    lower_area.x <- c(x$quantiles$x[y <= tel], max(y), tel)
    lower_area.y <- c(x$quantiles$q[y <= tel]^2, 1, 1)
  } else {
    lower_area.x <- x$quantiles$x[y <= tel]
    lower_area.y <- x$quantiles$q[y <= tel]^2
  }

  g <- ggplot2::ggplot(data.frame(x=x$quantiles$x, y=x$quantiles$q), mapping=ggplot2::aes(x=x, y=y)) +
    ggplot2::theme_light() +
    ggplot2::geom_line(size=1.3)
  if (any(y<=tel)) {
    g <- g + ggplot2::geom_ribbon(data=data.frame(x=lower_area.x, y=lower_area.y), mapping=ggplot2::aes(ymin=0, ymax=y), fill="gray90", col="gray60")
  }
    if (any(y>=tel)) {
    g <- g + ggplot2::geom_ribbon(data=data.frame(x=upper_area.x, y=upper_area.y), mapping=ggplot2::aes(ymin=1-y, ymax=1), fill="gray90", col="gray60")
  }
  g <- g + ggplot2::geom_line(data=data.frame(x=c(0, tel), y=c(0,0)), size=1.3, col='darkolivegreen4') +
    ggplot2::geom_line(data=data.frame(x=c(tel, tel), y=c(0,1)), size=1.3, col='darkolivegreen4') +
    ggplot2::geom_line(data=data.frame(x=c(tel, max(c(lower_area.x, upper_area.x))), y=c(1,1)), size=1.3, col='darkolivegreen4', linetype="solid") +
    ggplot2::xlab("Power [MW]") +
    ggplot2::ylab("Cumulative Distribution") +
    ggplot2::ggtitle(paste("CRPS:", round(CRPS(x, tel), 2), "MW"))

  plot(g)
}


#' Check probabilistic forecast class
is.prob_forecast <- function(x) inherits(x, "prob_forecast")

#' Register generic sample function
#' @param x A prob_forecast object
get_1d_samples <- function(x, ...) {
    UseMethod("get_1d_samples",x)
}

#' Register generic joint density function
#' @param x An n-dimension prob_forecast object
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
CVAR <- function(x, ...) {
  UseMethod("CVAR",x)
}


#' Register generic plot function
#' @param x A prob_forecast object
plot_pdf <- function(x, ...) {
  UseMethod("plot_pdf",x)
}

error_check_calc_quantiles_input <- function(quantiles){
  if (!(all(quantiles > 0 & quantiles < 1))) stop('Bad input. All quantiles must be in (0,1).')
}

# Methods for aggregate probabilistic forecast class using vine copulas
#------------------------------------------------------------------------------


#' Initialize a probabilistic power forecast for a specific time point, using an n-dimensional vine copula.
#' Assumes training data already captures differences in magnitude (i.e., power rating) amongst sites.
#' The "copula" field of the data.input list can be a matrix of training data [ntrain x nsites] OR a pre-trained vinecop model
#' The "marginals" field of the data.input list can be a matrix of training data [ntrain x nsites] OR a list of single-site prob_forecast objects
#'
#' @param data.input A list of "copula" and "marginals" inputs
#' @param location A string
#' @param time A lubridate time stamp
#' @param training_transform_type Transform of training data into uniform domain (see marg_transform "cdf.method")
#' @param results_transform_type Transform of copula results back into variable domain (see marg_transform "cdf.method")
#' @param n An integer, number of copula samples to take
#' @param samples.u (optional) A precalculated set of n-dimensional CDF samples from rvinecop
#' @param ... optional arguments to the marginal estimator
#' @return An n-dimensional probabilistic forecast object from vine copulas
fc_vine <- function(data.input, location, time,
                                  training_transform_type="empirical", results_transform_type='empirical', n=1e4,
                                  samples.u=NA, ...) {
  if (length(names(data.input))==0 | !all(names(data.input) == c("copula", "marginals"))) stop("data.input must be a list with 'copula' and 'marginals' inputs")
  if (!is.numeric(n)) stop('n (number of samples) must be an integer.')
  if (training_transform_type == "precalcbma" | results_transform_type == "precalcbma") {
    if (!all(sapply(data.input$marginals, is.prob_forecast))) stop("Single-site forecasts must be provided to use precalculated marginals. ")
    d <- length(data.input$marginals)
  } else {
    if (class(data.input$marginals)!='matrix') stop('Input data must be a matrix if marginals are to be calculated.')
    if (dim(data.input$marginals)[2] < 2) stop('Training data from more than 1 site required for vine copula forecast.')
    d <- dim(data.input$marginals)[2]
  }

  tr <- calc_transforms(data.input$marginals, d, training_transform_type, results_transform_type, ...)
  training_transforms <- tr$training
  results_transforms <- tr$results

  if ("vinecop" %in% class(data.input$copula)) {
    model <- data.input$copula
  } else {
    uniform_dat <- sapply(seq_along(training_transforms),
                          FUN=function(i) {to_uniform(training_transforms[[i]], data.input$copula[,i])}, simplify = "array")
    model <- rvinecopulib::vinecop(uniform_dat, family_set="all")
  }

  # Initialize probabilistic forecast
  dat <- list(training_transforms = training_transforms,
              results_transforms = results_transforms,
              location = location,
              time = time,
              model = model,
              d = d,
              n=n
              )
  x <- structure(dat, class = c("prob_forecast", "fc_vine"))

  # Complete probabilistic forecast by sampling and aggregating
  x$quantiles <- calc_quantiles(x, samples.u=samples.u)

  return(x)
}

#' Check class
is.fc_vine <- function(x) inherits(x, "fc_vine")

#' Calculate lists of variable-to-uniform domain transforms for all dimensions
#'
#' @param dat training data matrix
#' @param d Number of dimensions
#' @param training_transform_type Transform of training data into uniform domain (see marg_transform "cdf.method")
#' @param results_transform_type Transform of copula results back into variable domain (see marg_transform "cdf.method")
#' @param ... Optional arguments to marg_transform
#' @return list of "training" and "results" transforms to use.
calc_transforms <- function(dat, d, training_transform_type, results_transform_type, ...) {
  training <- lapply(seq_len(d), FUN=get_transform_with_unique_xmin_max, dat=dat, cdf.method=training_transform_type, ...)
  # Results transforms must be subsequently updated if desired
  if (results_transform_type==training_transform_type) {
    results <- training
  } else {
    results <- lapply(seq_len(d), FUN=get_transform_with_unique_xmin_max, dat=dat, cdf.method=results_transform_type, ...)
  }
  return(list('training'=training, 'results'=results))
}

#' Subfunction for calc_transforms to cycle thorugh xmin and xmax if they are given uniquely for each dimension
#'
#' @param idx Column index of dat
#' @param dat training data matrix over all the dimensions
#' @param cdf.method marg_transform cdf.method
#' @param ... Optional arguments to marg_transform, including potentially xmin or xmax in either scalar or vector form
get_transform_with_unique_xmin_max <- function(idx, dat, cdf.method, ...) {
  args <- list(...)
  # Use unique xmin/xmax values if vectors are given
  if ('xmin' %in% names(args) & length(args[['xmin']]) > 1) {args[['xmin']] <- args[['xmin']][idx]}
  if ('xmax' %in% names(args) & length(args[['xmax']]) > 1) {args[['xmax']] <- args[['xmax']][idx]}
  # Subset the proper column if data is a matrix or item if it is a list
  data <- if (class(dat)=='matrix') {dat[,idx]} else {dat[[idx]]}
  return(do.call(marg_transform, c(list(data, cdf.method), args))) # repackage arguments into single list for do.call
}

#' Sample the vine copula model and sum to calculate samples of the univariate, aggregate power forecast
#'
#' @return A column matrix of aggregate powers
get_1d_samples.fc_vine <- function(x, samples.u=NA) {
  # Get new samples from the copula, if needed.
  if (!(is.numeric(samples.u))) {
    samples.u <- rvinecopulib::rvinecop(x$n, x$model)
  }

  samples.xs <- sapply(seq_len(length(x)), FUN=function(i) return(from_uniform(x$results_transforms[[i]], samples.u[,i])), simplify="array")

  samples.x <- rowSums(samples.xs)
  return(samples.x)
}

#' Calculate forecast quantiles from samples of the vine copula
#'
#' @param x fc_vine object
#' @param samples (optional) previously obtained 1-D samples to use instead of new sampling, e.g. for coordination with cVaR calculation
#' @param samples.u (optional) previously obtained n-d samples to use instead of new sampling, e.g. for static vine copula model.
#' samples overrides samples.u if both are given
#' @param quantiles Sequence of quantiles in (0,1)
#' @return A list of q, the quantiles on [0, 1], and x, the estimated values
calc_quantiles.fc_vine <- function(x, samples=NA, samples.u=NA, quantiles=seq(0.001, 0.999, by=0.001)) {
  error_check_calc_quantiles_input(quantiles)

  if (!(is.numeric(samples))) {samples <- get_1d_samples(x, samples.u)}

  # Do I want this to be type 4 instead?
  xseq <- stats::quantile(samples, probs=quantiles, type=1, names=TRUE)

  return(list(x=xseq, q=quantiles))
}

#' Calculate VaR and CVaR from sampled data. CVaR calculation is done directly from the samples, rather than estimated from a fitted distribution.
#'
#' @param x fc_vine object
#' @param samples (optional) previously obtained samples to use instead of new sampling, e.g. for coordination with quantiles calculation
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return list of var, cvar
CVAR.fc_vine <- function(x, samples=NA, epsilon=c(0.05, 0.95)) {
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
#' @param x A fc_vine object
#' @param k Integer or vector of number of samples to take along each dimension. If given an integer, all dimensions have the same number of samples.
#' @return an joint probability density array
get_joint_density_grid.fc_vine <- function(x, k=100) {
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
#' @param x A fc_vine object
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
#' @param x fc_vine object
plot_pdf.fc_vine <- function(x, cvar=FALSE, epsilon=c(0.05, 0.95)) {
  # Assume data is power or irradiance and must be non-negative
  samples <- get_1d_samples(x)
  epdf <- stats::density(samples, from=0)
  plot(epdf, xlab='Power [MW]', ylab='Probability',
       main='Estimated probability density', sub = paste("Location: ", x$location, ", Time:", x$time))

  if (cvar){# Color in tails above/below desired epsilon's
    cvar_info <- CVAR(x, samples=samples, epsilon=epsilon)
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


# Methods for univariate probabilistic forecast class
#------------------------------------------------------------------------------

#' Quality control 1-dimension input data
#'
#' @param dat A numeric vector of data
#' @return data without NA's
qc_input <- function(dat) {
  if (!is.vector(dat)) stop("Input data must be a vector.")
  # Quality control for NA's
  dat <- dat[!is.na(dat)]
  if (length(dat) < 2) stop("Input data must have at least 2 non-NA values. ")
  return(dat)
}

#' Initialize a univariate probabilistic power forecast for a specific time point by ranking ensemble members.
#' Its rank quantiles are the basic quantiles estimated from the ensemble members; the quantiles are the rank quantiles iterpolated
#' to the quantiles of interest (i.e., 10%, 20%, etc...)
#'
#' @param data.input A numeric vector of ensemble members
#' @param location A string
#' @param time A lubridate time stamp
#' @param max_power Maximum power for normalizing forecast to [0,1]
#' @return A probabilistic forecast object
fc_binned <- function(data.input, location, time, max_power, ...) {
  members <- qc_input(data.input)

  # Initialize probabilistic forecast
  dat <- list(location = location,
              rank_quantiles = list(x=c(0, sort(members), max_power), y=(0:(length(members)+1))/(length(members)+1)),
              time = time,
              d = 1)

  x <- structure(dat, class = c("prob_forecast", "fc_binned"))
  x$quantiles <- calc_quantiles(x)

  return(x)
}

#' Check class
is.fc_binned <- function(x) inherits(x, "fc_binned")

#' Calculate forecast quantiles
#'
#' @param x fc_binned object
#' @param quantiles Sequence of quantiles in (0,1)
#' @return A list of q, the quantiles on [0, 1], and x, the estimated values
calc_quantiles.fc_binned <- function(x, quantiles=seq(0.001, 0.999, by=0.001)) {
  error_check_calc_quantiles_input(quantiles)
  xseq <- stats::approx(x=x$rank_quantiles$y,  y=x$rank_quantiles$x, xout=quantiles)$y

  return(list(x=xseq, q=quantiles))
}


#' Not implemented
#'
#' @param x fc_binned object
plot_pdf.fc_binned <- function(x) {
  stop("Not implemented for binned probability forecasts.")
}

# ---------------------------------------------------------------------------------------------

#' Initialize a forecast using an empirical (stepped, not interpolated) CDF of a training set.
#' Can be applied for a climatology, raw ensemble, or persistence ensemble type forecast.
#'
#' @param data.input A numeric vector of training data
#' @param location A string
#' @param time A lubridate time stamp
#' @return A probabilistic forecast object
fc_empirical <- function(data.input, location, time, ...) {
  data.input <- qc_input(data.input)

  # Initialize probabilistic forecast
  dat <- list(location = location,
              time = time,
              d = 1)

  x <- structure(dat, class = c("prob_forecast", "fc_empirical"))
  x$quantiles <- calc_quantiles(x, telemetry=data.input)

  return(x)
}

#' Check class
is.fc_empirical <- function(x) inherits(x, "fc_empirical")

#' Calculate forecast quantiles
#'
#' @param x fc_empirical object
#' @param telemetry A vector of telemetry to generate the empirical cdf from
#' @param quantiles Sequence of quantiles in (0,1)
#' @return A list of q, the quantiles on [0, 1], and x, the estimated values
calc_quantiles.fc_empirical <- function(x, telemetry=FALSE, quantiles=seq(0.001, 0.999, by=0.001)) {
  if (all(!telemetry)) stop("telemetry is a required input.")
  error_check_calc_quantiles_input(quantiles)
  xseq <- stats::quantile(telemetry, probs=quantiles, names=F, na.rm=T, type=1)

  return(list(x=xseq, q=quantiles))
}

#' Not implemented
#'
#' @param x fc_empirical object
plot_pdf.fc_empirical <- function(x) {
  stop("Not implemented for empirical forecasts.")
}

# ---------------------------------------------------------------------------------------------

#' Initialize a univariate probabilistic power forecast for a specific time point using Bayesian model averaging (BMA)
#'
#' @param data.input A vector of ensemble members
#' @param location A string
#' @param time A lubridate time stamp
#' @param model A pre-fit BMA model from bma_ens_models
#' @param max_power Maximum power for normalizing forecast to [0,1]
#' @param bma_distribution One of "beta", "truncnorm" --> determines the type of distribution used for each member's kernel dressing
#' @param ... Additional parameters
#' @return A probabilistic forecast object
fc_bma <- function(data.input, location, time, model, max_power, bma_distribution, ...) {
  # Sanity check inputs; skip if model is missing
  if (all(is.na(model))) return(NA)
  if (!(bma_distribution %in% c("beta", "truncnorm"))) stop("bma_distribution not recognized. Must be 'beta' or 'truncnorm'.")

  # Use specific quality control to handle missing members and adjust model weights accordingly
  model <- qc_bma_input(data.input, model)

  # Initialize probabilistic forecast
  dat <- list(location = location,
              time = time,
              d = 1,
              model=model,
              members=data.input,
              max_power=max_power,
              bma_distribution=bma_distribution
  )
  x <- structure(dat, class = c("prob_forecast", "fc_bma"))

  # Complete probabilistic forecast by sampling and aggregating
  dc_model <- get_discrete_continuous_model(x)
  x$geometry_codes <- dc_model$geometries
  x$quantiles <- calc_quantiles(x, model=dc_model)
  return(x)
}

#' Check class
is.fc_bma <- function(x) inherits(x, "fc_bma")

#' Specific quality control handling for BMA forecasts, for missing members that have already been trained in the model.
#' Reallocates member weights if necessary
qc_bma_input <- function(members, model) {
  if (length(members[!is.na(members)]) < 2) stop("Input data must have at least 2 non-NA values.")

  missing <- is.na(members)
  missing_weight <- sum(model$w[missing], na.rm=T)
  model$w[missing] <- 0
  model$w <- model$w/sum(model$w, na.rm=T)
  return(model)
}

#' Calculate forecast quantiles
#' @param x fc_bma object
#' @param quantiles Sequence of quantiles in (0,1)
#' @return A list of q, the quantiles on [0, 1], and x, the estimated values
calc_quantiles.fc_bma <- function(x, model=NA, quantiles=seq(0.001, 0.999, by=0.001)) {
  error_check_calc_quantiles_input(quantiles)

  if (all(is.na(model))) model <- get_discrete_continuous_model(x)

  xseq <- stats::approx(x=model$cdf, y=model$xseq, xout=quantiles, yleft=0)$y # yright=x$max_power

  d <- stats::approx(x=model$xseq, y=model$pdf, xout=xseq, yleft=0, yright = 0)$y

  return(list(x=xseq, q=quantiles, d=d))
}

#' @param x fc_bma object
#' @param xseq Vector of x values in [0,1] to evaluate
#' @param discrete Boolean of whether to return density with a discrete component (approximate) or with the continuous estimation (congruent with CDF)
#' @return A list of the discrete components (PoC) and continuous density (pdf) and distribution (cdf) for each member
get_discrete_continuous_model <- function(x, xseq=seq(0, 1, 0.0001), discrete=F) {
  i.thresh <- max(which(xseq < x$model$percent_clipping_threshold))

  shape_params <- get_shape_params(x)

  # Get geometry codes
  codes <- mapply(FUN=get_beta_distribution_geometry_code, x$bma_distribution, shape_params$param1s, shape_params$param2s)
  geometries <- list("U type"=sum(codes==1), "Reverse J"=sum(codes==2), "J-type"=sum(codes==3), "Upside-down U"=sum(codes==4), "Missing"=sum(codes==0))


  # Get sequence of beta cumulative distribution, by ensemble member. Returns matrix of [xseq x members]
  cdf_member <- mapply(cdf_subfunction, shape_params$param1s, shape_params$param2s, shape_params$PoC, x$model$w,
                      MoreArgs = list(xseq=xseq, i.thresh=i.thresh, bma_distribution=x$bma_distribution, max_power=x$max_power))

  # Get sequence of beta probability densities, by ensemble member. Returns matrix of [xseq x members]
  pdf_member <- mapply(pdf_subfunction, shape_params$param1s, shape_params$param2s, shape_params$PoC, x$model$w,
                      MoreArgs = list(xseq=xseq, i.thresh=i.thresh, discrete=discrete, bma_distribution=x$bma_distribution, max_power=x$max_power))

  # Calculate overall distribution
  PoC_total <- sum(x$model$w*shape_params$PoC, na.rm=T)
  pdf_total <- apply(X=pdf_member, MARGIN=1, FUN=function(row) sum(row, na.rm=T))
  cdf_total <- apply(X=cdf_member, MARGIN=1, FUN=function(row) sum(row, na.rm=T))
  # Ensure pbeta sums to 1 after weighting and summing
  if (max(cdf_total) > 1) cdf_total <- cdf_total/max(cdf_total)

  # Invert normalization of beta components back to [0, max power]
  return(list(PoC=PoC_total, pdf=pdf_total/x$max_power, cdf=cdf_total, xseq=xseq*x$max_power, geometries=geometries,
              members=list(PoC=shape_params$PoC, pdf=pdf_member/x$max_power, cdf=cdf_member, codes=codes)))
}

# Get cumulative distribution for individual ensemble member
#' @param param1 First parameter: alpha for a beta distribution; mean for a truncnorm
#' @param param2 First parameter: beta for a beta distribution; sd for a truncnorm
#' @param bma_distribution "beta" or "truncnorm"
cdf_subfunction <- function(param1, param2, poc, w, xseq, i.thresh, bma_distribution, max_power) {
  if (bma_distribution=="beta") {
    p <- stats::pbeta(xseq[1:i.thresh], param1, param2)
  } else if (bma_distribution=="truncnorm") {
    p <- truncnorm::ptruncnorm(xseq[1:i.thresh], a=0, b=max_power, mean=param1, sd=param2)
  } else stop(paste("bma_distribution not recognized. Given", bma_distribution))

  distribution.continuous <- (1-poc)*p/max(p)
  dx <- xseq[length(xseq)]-xseq[i.thresh]

  # Model discrete component as a linear increase over the clipping bandwidth: y = mx + b, where b=1-poc, m=poc/clipping bandwidth
  distribution.discrete <- (1-poc) + seq_along(xseq[(i.thresh+1):length(xseq)])*diff(xseq)[1]*poc/dx
  return(w * c(distribution.continuous, distribution.discrete))
}

#' Get probability density for individual ensemble member
#' Includes continuous estimation of discrete component, starting at the threshold and extending to 1
#' @param bma_distribution "beta" or "truncnorm"
#' @param param1 First parameter: alpha for a beta distribution; mean for a truncnorm
#' @param param2 First parameter: beta for a beta distribution; sd for a truncnorm
pdf_subfunction <- function(param1, param2, poc, w, xseq, i.thresh, discrete, bma_distribution, max_power) {
  if (discrete) {
    if (bma_distribution=="beta") {
      density <- stats::dbeta(xseq, param1, param2)
    } else if (bma_distribution=="truncnorm") {
      density <- truncnorm::dtruncnorm(xseq, a=0, b=max_power, mean=param1, sd=param2)
    } else stop(paste("bma_distribution not recognized. Given", bma_distribution))
    return((1-poc)*w*density)
  } else {

    if (bma_distribution=="beta") {
      distribution_max <- stats::pbeta(xseq[i.thresh], param1, param2)
      density <- stats::dbeta(xseq[1:i.thresh], param1, param2)
    } else if (bma_distribution=="truncnorm") {
      distribution_max <- truncnorm::ptruncnorm(xseq[i.thresh], a=0, b=max_power, mean=param1, sd=param2)
      density <- truncnorm::dtruncnorm(xseq[1:i.thresh], a=0, b=max_power, mean=param1, sd=param2)
    } else stop(paste("bma_distribution not recognized. Given", bma_distribution))

    density.continuous <- (1-poc)*density/distribution_max
    dx <- xseq[length(xseq)]-xseq[i.thresh]

    # Model discrete component as a linear increase over the clipping bandwidth: y = mx + b, where b=1-poc, m=poc/clipping bandwidth
    density.discrete <- poc/(xseq[length(xseq)]-xseq[i.thresh]) * rep(1, times=length(xseq[(i.thresh+1):length(xseq)]))
    return(w * c(density.continuous, density.discrete))
  }
}

#' Get shape parameters for either a beta or a truncated normal distribution
#' @param x A fc_bma object
#' @return a list of param1s (alpha for a beta distribution; mean for a truncnorm), param2s (beta for a beta distribution; sd for a truncnorm), and POC (probability of clipping)
get_shape_params <- function(x) {
   # Normalize to [0,1]
   members.norm <- sapply(x$members, function(m) ifelse(m <= x$max_power, m/x$max_power, 1))

   # Get parameters for individual ensemble members
   PoC <- mapply(get_poc, members.norm, x$model$A0, x$model$A1, x$model$A2, MoreArgs=list(A_transform=x$model$A_transform))
   rhos <- mapply(get_rho, members.norm, x$model$B0, x$model$B1, MoreArgs = list(B_transform=x$model$B_transform))
   sigmas <- mapply(get_sigma, rhos, MoreArgs = list(C0=x$model$C0))

   if (x$bma_distribution=="beta") {
     gammas <- mapply(get_gamma, rhos, sigmas)

     # Truncate gammas if needed at a J or reverse-J shape to avoid U-shaped distributions
     # variance is inversely proportional to gamma, so gamma min leads to variance max
     gamma_min <- sapply(rhos, function(r) ifelse(r<=0.5, 1/(1-r), 1/r))
     gammas[gammas < gamma_min & !is.na(gammas)] <- gamma_min[gammas < gamma_min & !is.na(gammas)]

     param1s <- rhos * gammas # alphas
     param2s <- gammas * (1-rhos) #betas
   } else if (x$bma_distribution=="truncnorm") {
     param1s <- rhos
     param2s <- sigmas
   } else stop(paste("bma_distribution not recognized. Given", x$bma_distribution))

   return(list(param1s=param1s, param2s=param2s, PoC=PoC))
 }

# Broadest geometry categories. Most important is identifying and eliminating U-shaped.
#' 0 -> Missing
#' 1 -> U-type
#' 2 -> Reverse J-type
#' 3 -> J-type
#' 4 -> Upside-down U type
get_beta_distribution_geometry_code <- function(bma_distribution, alpha, beta) {
  if (bma_distribution == "truncnorm") {return(0)} # Not applicable

  if (is.na(alpha) | is.na(beta)) {return(0)}
  if (alpha < 1 & beta < 1) {return(1)}
  else if (alpha < 1 & beta >= 1) {return(2)}
  else if (alpha >= 1 & beta < 1) {return(3)}
  else if (alpha >= 1 & beta >= 1) {return(4)}
  else stop(paste("uncoded combination: alpha=", alpha, ", beta=", beta, sep=''))
}

#' Plot APPROXIMATION of the BMA probability density function, including the member component contributributions
#' @param discrete Boolean on whether to plot a discrete component as approximation or the continuous equivalent above the clipping threshold
plot_pdf.fc_bma <- function(x, actual=NA, ymax=NA, normalize=F, discrete=T) {

  model <- get_discrete_continuous_model(x, discrete=discrete)

  # Avoid Inf's on the boundaries
  xrange <- 2:(length(model$xseq)-1)
  ymax <- ifelse(is.na(ymax), max(c(max(model$pdf[xrange])*ifelse(normalize, x$max_power, 1), model$PoC))*1.1, ymax)

  g <- ggplot2::ggplot(data.frame(x=model$xseq[xrange]/ifelse(normalize, x$max_power, 1),
                                  y=model$pdf[xrange]*ifelse(normalize, x$max_power, 1)),
                       mapping=ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_line(size=1.3) +
    ggplot2::xlab(ifelse(normalize, "Normalized Power [MW]", "Power [MW]")) +
    ggplot2::ylab("Probability Density") +
    ggplot2::geom_point(data=data.frame(x=x$members/ifelse(normalize, x$max_power, 1), y=ymax), col="black", fill="grey", alpha=0.5, shape=21, size=3)
  if (discrete) {
    # Lollipop style segment for discrete component
    g <- g + ggplot2::geom_line(data=data.frame(x=c(ifelse(normalize, 1, x$max_power), ifelse(normalize, 1, x$max_power)),
                                       y=c(0,model$PoC)), size=1.3) +
      ggplot2::geom_point(data=data.frame(x=c(ifelse(normalize, 1, x$max_power)), y=c(model$PoC)), size=4, col="black", fill="black", shape=21)
  }

  if (!is.na(actual)) {
    g <- g + ggplot2::geom_line(data=data.frame(x=c(actual/ifelse(normalize, x$max_power, 1), actual/ifelse(normalize, x$max_power, 1)),
                                                y=c(0, ymax)), linetype="dashed")
  }

  for (i in seq_len(dim(model$members$pdf)[2])) {
    g <- g + ggplot2::geom_line(data.frame(x=model$xseq[xrange]/ifelse(normalize, x$max_power, 1),
                                           y=model$members$pdf[xrange,i]*ifelse(normalize, x$max_power, 1)),
                                mapping=ggplot2::aes(x=x, y=y), col="black")
  }

  plot(g)
}



# ---------------------------------------------------------------------------------------------

#' Initialize a univariate probabilistic power forecast for a specific time point using kernel density estimation. See kde_methods.R for more details
#'
#' @param data.input A vector of ensemble members
#' @param location A string
#' @param time A lubridate time stamp
#' @param cdf.method KDE method selection, see kde_methods.R for details
#' @param ... Additional parameters passed on the KDE method
#' @return A probabilistic forecast object
fc_kde <- function(data.input, location, time, cdf.method='geenens', ...) {
  members <- qc_input(data.input)

  func <- marginal_lookup(cdf.method)
  # Get selected KDE
  model <- func(members, ...)

  # Initialize probabilistic forecast
  dat <- list(location = location,
              time = time,
              d = 1,
              model=model
              )
  x <- structure(dat, class = c("prob_forecast", "fc_kde"))

  # Complete probabilistic forecast by sampling and aggregating
  x$quantiles <- calc_quantiles(x)
  return(x)
}

#' Check class
is.fc_kde <- function(x) inherits(x, "fc_kde")

#' Calculate forecast quantiles from KDE
#'
#' @param x fc_kde object
#' @param quantiles Sequence of quantiles in (0,1)
#' @return A list of q, the quantiles on [0, 1], and x, the estimated values
calc_quantiles.fc_kde <- function(x, quantiles=seq(0.001, 0.999, by=0.001)) {
  error_check_calc_quantiles_input(quantiles)
  xseq <- stats::approx(x=x$model$u,  y=x$model$x, xout=quantiles)$y

  return(list(x=xseq, q=quantiles))
}

#' Calculate VaR and CVaR by trapezoidal integration.
#'
#' @param x fc_kde object
#' @param epsilon Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)
#' @return list of var, cvar
CVAR.fc_kde <- function(x, epsilon=c(0.05, 0.95)) {
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
#' @param x fc_kde object
plot_pdf.fc_kde <- function(x) {
  plot(x$model$x, x$model$d, xlab='Power [MW]', ylab='Probability density', sub = paste("Location: ", x$location, ", Time:", x$time),
       type='l', lwd=2)
}
