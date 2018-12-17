#' Initialize a time series of power forecasts. Assumes training data already captures
#' differences in magnitude (i.e., power rating) amongst sites. Forecast is NA for times when sun is down.
#'
#' @param x A matrix or array of training data of dimensions [time x ntrain x nsites] for n-dimensional forecast or [time x ntrain] for a 1-dimensional forecast
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'site', 'region', 'total'
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine' (irrelevant if scale == 'site')
#' @param ... optional arguments to the prob_forecast object
ts_forecast <- function(x, start_time, time_step, scale, location, method, ...) {
  # Check inputs
  if (!is.array(x)) stop("Bad input. x must be an array.")
  if (!((length(dim(x))==2 & tolower(scale) %in% c("site", "s")) | (length(dim(x))==3 & tolower(scale) %in% c("region", "total", "r", "t")))){
    stop("Data and scale mis-match. x must be a 2-dimensional array for 1-dimensional forecasts or a 3-dimensional array for n-dimensional forecasts.")
  }

  sun_up <- apply(x, MARGIN=1, FUN=check_sunup)
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
  class_type <- get_forecast_class(scale, method)
  sub_func <- function(i, time, sun_up, x, location, ...){
    if (sun_up) {
      if (class_type$dim == '1') {
        return(class_type$class(x[i,], location=location, time=time, ...))
      } else {return(class_type$class(x[i,,], location=location, time=time, ...))}
    } else {return(NA)}}
  return(mapply(sub_func, seq_len(dim(x)[1]), start_time + (seq_len(dim(x)[1])-1)*lubridate::dhours(time_step), sun_up,
                MoreArgs=list(x=x, location=location, ...=...), SIMPLIFY=FALSE))
}

#' Look-up function of the forecast class type
#'
#' @param scale One of 'site', 'region', 'total'
#' @param method One of 'rank' if scale is 'site', else one of gaussian', 'empirical', 'vine'
#' @return A function to initialize a forecast of the desired type
get_forecast_class <- function(scale, method){
  if (tolower(scale) %in% c("site", "s")){
    cls <- switch(tolower(method),
                  "rank" = prob_1d_rank_forecast,
                  "r" = prob_1d_rank_forecast,
                  "kde" = prob_1d_kde_forecast,
                  "k" = prob_1d_kde_forecast,
                  stop(paste('Forecast type', method, 'not recognized for single-site forecasts.', sep=' ')))
    d <- '1'
  } else if (tolower(scale) %in% c("region", "total", "r", "t")){
    cls <- switch(tolower(method),
           "vine" = prob_nd_vine_forecast,
           "gaussian" = prob_nd_gaussian_forecast,
           "empirical" = prob_nd_empirical_forecast,
           "v" = prob_nd_vine_forecast,
           "g" = prob_nd_gaussian_forecast,
           "e" = prob_nd_empirical_forecast,
           stop(paste('Forecast type', method, 'not recognized for multi-dimensional forecasts.', sep=' ')))
    d <- 'n'
  } else stop(paste('Forecast scale', scale, 'not recognized.', sep=' '))
  return(list(class=cls, dim=d))
  }

#' Calculate number of steps in time series forecast
length.ts_forecast <- function(x) {
  return(length(x$forecasts))
}

#' Plot fan plot of the time-series probabilistic forecast
#'
#' @param x A ts_forecast object
#' @param tel vector of telemetry (optional)
plot.ts_forecast <- function(x, ..., tel=NA) {
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
  if (any(!is.na(tel))) {
    actuals_time_step <- x$time_step*length(x)/length(tel)
    graphics::lines(seq(from=actuals_time_step, length.out=length(tel), by=actuals_time_step), tel, col='chocolate3', lwd=2)
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

# --------------------------------------------------------------------------------------------
# Forecast evaluation methods

#' Preprocess for metrics evaluations, for when telemetry is at finer time resolution than the forecast.
#' The forecast and telemetry can be at different time resolutions, so long as telemetry is a multiple of the forecast.
#'
#' @param tel A vector of the telemetry values
#' @param fc A vector of data from the time series forecast
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @return list of the telemetry and forecast data vectors, of equal length
equalize_telemetry_forecast_length <- function(tel, fc, agg=TRUE) {
  # Pre-process vectors to same lengths
  if (length(tel) != length(fc) & length(tel) %% length(fc) > 0) stop("Telemetry length must be equal to or a multiple of forecast length.")
  tel_2_fc <- length(tel)/length(fc)
  if (agg & tel_2_fc > 1) {  # Aggregate telemetry
    tel <- sapply(seq_along(fc), function(i) {return(sum(tel[(tel_2_fc*(i-1)+1):(tel_2_fc*i)])/tel_2_fc)})
  } else {  # Expand forecast data
    fc <- rep(fc, each=tel_2_fc)
  }
  return(list(tel=tel, fc=fc, tel_2_fc=tel_2_fc))
}

#' Get statistics on the prevalence of NaN's in telemetry and sun-up/sun-down discrepancies
#'
#' @param tel A vector of the telemetry values
#' @param len_ts Number of timesteps in the time series forecast
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @param return A list of summary statistics
get_sundown_and_NaN_stats <- function(ts, tel, agg=TRUE) {
  res <- equalize_telemetry_forecast_length(tel, ts$sun_up, agg=agg)
  sun_up <- res$fc
  tel <- res$tel

  fc_sundown <- length(sun_up[sun_up==FALSE])
  fc_sunup <- length(sun_up[sun_up==TRUE])
  tel_sunup_NaN <- length(tel[is.nan(tel) & sun_up==TRUE])
  tel_sunup <- length(tel[!is.nan(tel) & sun_up==TRUE])
  tel_sundown_0 <- length(tel[!is.nan(tel) & tel==0 & sun_up==FALSE])
  tel_sundown <- length(tel[!is.nan(tel) & tel!=0 & sun_up==FALSE])
  tel_sundown_NaN <- length(tel[is.nan(tel) & sun_up==FALSE])
  return(list("Sun-up forecasts"=fc_sunup,
              "Sun-down forecasts"=fc_sundown,
              "Telemetry missing when sun up"=tel_sunup_NaN,
              "Telemetry available when sun up"=tel_sunup,
              "Telemetry is 0 when sun forecasted down"=tel_sundown_0,
              "Telemetry is non-zero when sun forecasted down"=tel_sundown,
              "Telemetry is NaN when sun forecasted down"=tel_sundown_NaN
         ))
}

#' Get the average estimated CRPS (continuous ranked probability score) for the forecast;
#' the score at each time point is estimated from sampled data.
#' CRPS characterizes calibration and sharpness together
#'
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
eval_avg_crps <-function(ts, tel){
  tel_e <- aggregate_telemetry(tel, length(ts))
  crps <- mapply(function(forecast, value){
    if (is.prob_forecast(forecast)){return(scoringRules::crps_sample(value, get_1d_samples(forecast)))}
    else {return(NA)} }, ts$forecasts, tel_e)
  return(list(sd=stats::sd(crps[ts$sun_up]), mean= mean(crps[ts$sun_up])))
}

#' Get Brier score at a certain probability of exceedance
#' This is calculated with a constant probability threshold, rather than a constant value threshold
#' I THINK THIS NEEDS TO BE EXPLORED, BECAUSE I DON'T KNOW THE EFFECT OF THE CONSTANT PROBABILITY THRESHOLD.
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param alpha Threshold probability of exceedance, numeric [0,1]
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @return the Brier score
eval_brier <- function(ts, tel, alpha, agg=TRUE) {
  if (alpha < 0 | alpha > 1) stop(paste("alpha must be [0,1], given ", alpha, '.', sep=''))
  thresholds <- get_quantile_time_series(ts, 100*(1-alpha))

  thresholds <- equalize_telemetry_forecast_length(tel, thresholds, agg=agg)$fc
  res <- equalize_telemetry_forecast_length(tel, ts$sun_up, agg=agg)
  sun_up <- res$fc
  tel <- res$tel

  indicator <- as.integer(tel[sun_up] >= thresholds[sun_up])
  return(sum(((1-alpha)-indicator)^2, na.rm = TRUE))
}

#' Get mean absolute error between the forecast median and the actual value
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @return the MAE value
eval_mae <-function(ts, tel, agg=TRUE) {
  medians <- get_quantile_time_series(ts, 50)
  sun_up <- equalize_telemetry_forecast_length(tel, ts$sun_up, agg=agg)$fc
  res <- equalize_telemetry_forecast_length(tel, medians, agg=agg)
  medians <- res$fc
  tel <- res$tel
  return(mean(abs(medians[sun_up]-tel[sun_up]), na.rm = TRUE))
}

#' Get average interval score, for an interval from alpha/2 to 1-alpha/2. Negatively oriented (smaller is better)
#' Characterizes sharpness, with a penalty for reliability
#' NEED TO EXPLORE THIS. UNSURE OF VALUE OF DOING SIMPLE AVERAGE OVER A HETEROSCEDASTIC PROCESS TO CHARACTERIZE SHARPNESS.
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @return the average IS value
eval_avg_is <-function(ts, tel, alpha, agg=TRUE) {
  x <- equalize_telemetry_forecast_length(tel, ts$sun_up, agg=agg)
  sun_up <- x$fc
  return(mean(sapply(which(sun_up), function(i) {if (agg) {j <- i} else {j <- floor((i-1)/x$tel_2_fc)+1}
    calc_is(ts$forecasts[[j]], x$tel[i], alpha=alpha)}), na.rm=TRUE))
}

#' Plot diagonal line diagram of quantiles + observations
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
plot_reliability <- function(ts, tel, agg=TRUE) {
  x <- get_quantile_reliability(ts, tel, agg=agg)

  graphics::plot(x$quantiles, x$quantiles, type="l", lty=2, xlab="Nominal",
                 ylab="Observed", main="To do: add uncertainty bars")
  graphics::lines(x$quantiles, cumsum(x$hit_rate), type='b', lty=1, pch=1)
}

#' Evaluate forecast reliability by evaluating the actual hit rate of the quantile bins
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @return list of the quantiles and their hit rates
get_quantile_reliability <- function(ts, tel, agg=TRUE) {
  x <- equalize_telemetry_forecast_length(tel, ts$sun_up, agg=agg)
  sun_up <- x$fc

  # Get the list of quantiles that have been evaluated, and add top limit at 1
  quants <- c(sapply(names(ts$forecasts[[min(which(sapply(ts$forecasts, FUN=is.prob_forecast)))]]$quantiles),
         function(x) {as.numeric(gsub("%", "", x))/100}, USE.NAMES=FALSE), 1)

  # Find time-points where telemetry and forecast data is available
  indices <- which(!is.nan(x$tel) & sun_up==TRUE)
  q_idx_subfunc <- function(i) {
    if (agg) {j <- i} else {j <- floor((i-1)/x$tel_2_fc)+1}
    list_idx <- which(ts$forecasts[[j]]$quantiles > x$tel[i]) # List of quantile indices above telemetry value
    if (length(list_idx) > 0) {idx <- min(list_idx)} else{idx <- length(quants)} # Pick lowest quantile, or 100th percentile if it falls outside distribution
    return(idx)
  }
  q_indices <- sapply(indices, q_idx_subfunc) # Get the indices of the corresponding quantile for each time points
  counts <- sapply(seq_along(quants), function(i) return(sum(q_indices == i)), USE.NAMES = FALSE)

  hit_rate <- counts/length(indices)
  return(list(quantiles=quants, hit_rate=hit_rate))
}
