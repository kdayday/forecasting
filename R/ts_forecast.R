#' Initialize a time series of power forecasts. Assumes training data already captures
#' differences in magnitude (i.e., power rating) amongst sites. Forecast is NA for times when sun is down.
#'
#' @param x A list of training data the length of the time-series. Each element should be [ntrain x 1] matrix of training data (for scale=='site')
#' or a [ntrain x nsites] matrix for scale=='region' or 'total'
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'site', 'region', 'total'
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine' (irrelevant if scale == 'site')
#' @param ... optional arguments to the prob_forecast object
ts_forecast <- function(x, start_time, time_step, scale, location, method, ...) {
  # Check inputs
  if (!((is.vector(x[[1]]) & tolower(scale) %in% c("site", "s")) | (is.matrix(x[[1]]) & tolower(scale) %in% c("region", "total", "r", "t")))){
      stop("Data and scale mis-match. x should be a list the length of the time-series. Each element should be a vector of training data (for scale=='site') or a [ntrain x nsites] matrix for scale=='region' or 'total")
  }

  sun_up <- unlist(lapply(x, check_sunup))
  forecasts <- calc_forecasts(x, sun_up, start_time=start_time, time_step=time_step, scale=scale, location=location,
                              method=method, ...)

  dat <- list(start_time = start_time,
              scale = scale,
              location = location,
              time_step = time_step,
              forecasts = forecasts,
              sun_up = sun_up
  )
  structure(dat, class = "ts_forecast")
}

#' Check GHI/power matrix for any positive numbers
#'
#' @param x A matrix
#' @return A boolean
check_sunup <- function(x){
  return(any(x>0, na.rm=TRUE))
}

#' Calculate a time series list of power forecasts.
#'
#' @param x A list of training data the length of the time-series. Each element should be [ntrain x 1] matrix of training data (for scale=='site')
#' or a [ntrain x nsites] matrix for scale=='region' or 'total'
#' @param sun_up Logical vector the length of the time-series
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'site', 'region', 'total'
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine'
#' @param ... optional arguments to the prob_forecast, including marginal estimator arguments
#' @return A list of forecasts. Forecast is NA for times when sun is down.
calc_forecasts <- function(x, sun_up, start_time, time_step, scale, location, method, ...) {
  forecast_class <- get_forecast_class(scale, method)
  forecasts <- mapply(function(d, time, sun_up, ...){if (sun_up) {return(forecast_class(d[!is.na(d)], location=location, time=time, ...))} else {return(NA)}},
                      x, start_time + (seq_along(x)-1)*lubridate::dhours(time_step), sun_up, MoreArgs=list(...=...), SIMPLIFY=FALSE)
  return(forecasts)
}

#' Look-up function of the forecast class type
#'
#' @param scale One of 'site', 'region', 'total'
#' @param method One of 'rank' if scale is 'site', else one of gaussian', 'empirical', 'vine'
#' @return A function to initialize a forecast of the desired type
get_forecast_class <- function(scale, method){
  if (tolower(scale) %in% c("site", "s")){
    return(switch(tolower(method),
                  "rank" = prob_1d_rank_forecast,
                  "r" = prob_1d_rank_forecast,
                  "kde" = prob_1d_kde_forecast,
                  "k" = prob_1d_kde_forecast,
                  stop(paste('Forecast type', method, 'not recognized for single-site forecasts.', sep=' '))))
  } else if (tolower(scale) %in% c("region", "total", "r", "t")){
    return(switch(tolower(method),
           "vine" = prob_nd_vine_forecast,
           "gaussian" = prob_nd_gaussian_forecast,
           "empirical" = prob_nd_empirical_forecast,
           "v" = prob_nd_vine_forecast,
           "g" = prob_nd_gaussian_forecast,
           "e" = prob_nd_empirical_forecast,
           stop(paste('Forecast type', method, 'not recognized for multi-dimensional forecasts.', sep=' '))))
  } else stop(paste('Forecast scale', scale, 'not recognized.', sep=' '))
}

#' Calculate number of steps in time series forecast
length.ts_forecast <- function(x) {
  return(length(x$forecasts))
}

#' Plot fan plot of the time-series probabilistic forecast
#'
#' @param x A ts_forecast object
#' @param actuals list of realized values (optional)
plot.ts_forecast <- function(x, ..., actuals=NA) {
  # Hard-coding this for now -- will need to change if other percentiles are desired
  probs <- seq(from=0, to=1, length.out = length(ts$forecasts[[min(which(sapply(ts$forecasts, FUN=is.prob_forecast)))]]$quantiles))
  plotdata <- matrix(ncol=length(x), nrow=length(probs))
  for (i in seq_along(x$forecasts)) {
    if (is.prob_forecast(x$forecasts[[i]])) {
      plotdata[,i] <- unname(x$forecasts[[i]]$quantiles)
    } else plotdata[,i] <- 0
  }
  graphics::plot(NULL, xlim=c(0, length(x)*x$time_step), ylim=c(0, max(plotdata)), xlab="Time [Hrs]", ylab="Aggregate power [W]")
  fanplot::fan(plotdata, data.type='values', probs=probs, fan.col=colorspace::sequential_hcl,
               rlab=NULL, start=x$time_step, frequency=1/x$time_step)
  if (all(!is.na(actuals))) {
    graphics::lines(seq(from=x$time_step, length.out=length(x), by=x$time_step), actuals, col='chocolate3', lwd=2)
  }
}

#' Get a time-series of the forecast value at a particular quantile
#' @param x A ts_forecast object
#' @param alpha Quantile of interest, numeric [0, 100]
#' @return Numeric vector
get_quantile_time_series <- function(x, alpha) {
  q <- paste(alpha, '%', sep='')
  timeseries <- sapply(x$forecasts, FUN=function(forecast, q) {if (is.prob_forecast(forecast)){return(unname(forecast$quantiles[q]))} else return(0)}, q=q)
  if (any(is.na(timeseries[x$sun_up]))) stop(paste(q,'quantile not available'))
  return(timeseries)
}

#' Plot a time-series of the upper and lower tail CVAR's
#' @param x A ts_forecast object
plot_cvar_over_time <- function(x) {
  y1 <- sapply(x$forecasts, FUN=function(fc) {if (is.prob_forecast(fc)) {return(fc$cvar$low)} else {return(0)}})
  y2 <- sapply(x$forecasts, FUN=function(fc) {if (is.prob_forecast(fc)) {return(fc$cvar$high)} else {return(0)}})
  colors <- c("darkseagreen4", "darkorange3")
  graphics::plot((1:length(x))*x$time_step, y2, type = 'b', col=colors[1], xlab='Time', ylab='CVaR [W]')
  graphics::lines((1:length(x))*x$time_step, y1, type = 'b', col=colors[2])
  graphics::legend(x="topleft", legend=c("Upper tail", "Lower tail"), bty='n', col=colors, lty=1, lwd=2, y.intersp=2)
}

#' Get the average estimated CRPS (continuous ranked probability score) for the forecast;
#' the score at each time point is estimated from sampled data.
#' CRPS characterizes calibration and sharpness together
#'
#' @param ts A ts_forecast object
#' @param actuals A list of the realized values
eval_avg_crps <-function(ts, actuals){
  if (length(actuals) != length(ts)) {stop('Time-series forecast and actual values must have the same length.')}
  crps <- mapply(function(forecast, value){
    if (is.prob_forecast(forecast)){return(scoringRules::crps_sample(value, get_1d_samples(forecast)))}
    else {return(NA)} }, ts$forecasts, actuals)
  return(list(sd=stats::sd(crps[ts$sun_up]), mean= mean(crps[ts$sun_up])))
}

#' Get Brier score at a certain probability of exceedance
#' This is calculated with a constant probability threshold, rather than a constant value threshold
#' I THINK THIS NEEDS TO BE EXPLORED, BECAUSE I DON'T KNOW THE EFFECT OF THE CONSTANT PROBABILITY THRESHOLD.
#'
#' @param ts A ts_forecast object
#' @param actuals A list of the realized values
#' @param alpha Threshold probability of exceedance, numeric [0,100]
#' @return the Brier score
eval_brier <- function(ts, actuals, alpha) {
  if (alpha < 0 | alpha > 100) stop(paste("alpha must be [0,1], given ", alpha, '.', sep=''))
  thresholds <- get_quantile_time_series(ts, 100-alpha)
  indicator <- as.integer(actuals[ts$sun_up] >= thresholds[ts$sun_up])
  return(sum(((1-alpha/100)-indicator)^2))
}

#' Get mean absolute error between the forecast median and the actual value
#'
#' @param ts A ts_forecast object
#' @param actuals A list of the realized values
#' @return the MAE value
eval_mae <-function(ts, actuals) {
  medians <- get_quantile_time_series(ts, 50)
  return(mean(abs(medians[ts$sun_up]-actuals[ts$sun_up])))
}

#' Get average interval score, for an interval from alpha/2 to 1-alpha/2. Negatively oriented (smaller is better)
#' Characterizes sharpness, with a penalty for reliability
#' NEED TO EXPLORE THIS. UNSURE OF VALUE OF DOING SIMPLE AVERAGE OVER A HETEROSCEDASTIC PROCESS TO CHARACTERIZE SHARPNESS.
#'
#' @param ts A ts_forecast object
#' @param actuals A list of the realized values
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
#' @return the average IS value
eval_avg_is <-function(ts, actuals, alpha) {
  return(mean(mapply(calc_is, ts$forecasts[ts$sun_up], actuals[ts$sun_up], alpha=alpha)))
}

#' Plot diagonal line diagram of quantiles + observations
plot_reliability <- function(ts, actuals) {
  quants <- vector('numeric', length=length(ts$forecasts[[min(which(sapply(ts$forecasts, FUN=is.prob_forecast)))]]$quantiles))
  for (i in seq_along(ts)){
    if (is.prob_forecast(ts$forecasts[[i]])){
      indx <- min(which(ts$forecasts[[i]]$quantiles > actuals[i]))
      quants[indx] <- quants[indx] + 1
    }
  }
  quants <- quants/length(ts$forecasts[ts$sun_up])
  graphics::plot(seq(0, 1, along.with=quants), seq(0, 1, along.with=quants), type="l", lty=2, xlab="Nominal",
                 ylab="Observed", title="To do: add uncertainty bars")
  graphics::lines(seq(0, 1, along.with=quants), cumsum(quants), type='b', lty=1, pch=1)
}
