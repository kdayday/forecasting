#' Register ts_forecast constructor to allow for overloading
#' @param x A ts_forecast object
ts_forecast <- function(x, ...) {
  UseMethod("ts_forecast",x)
}

#' Construct a time series of power forecasts, using input training data. Assumes training data already captures
#' differences in magnitude (i.e., power rating) amongst sites. Forecast is NA for times when sun is down.
#'
#' @param x An array of training data of dimensions [time x ntrain x nsites] for n-dimensional forecast
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'region', 'total'
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine' (irrelevant if scale == 'site')
#' @param sun_up_threshold An absolute [MW] threshold on the ensemble members to remove dubious sunrise/sunset valud
#' @param MoreTSArgs An optional dictionary of time-series arguments to the forecast calculation
#' @param ... optional arguments to the prob_forecast object
ts_forecast.array <- function(x, start_time, time_step, scale, location, method, sun_up_threshold=0.5, MoreTSArgs=NA, ...) {
  # Check inputs
  if (!((length(dim(x))==3 & tolower(scale) %in% c("region", "total")))){
    stop("Data and scale mis-match. x must be a 3-dimensional array for n-dimensional forecasts.")
  }
  if (any(names(MoreTSArgs) %in% c("location", "time"))) stop("MoreTSArgs may not include location or time, which are reserved names.")

  new_ts_forecast(x, start_time, time_step, scale, location, method, sun_up_threshold, MoreTSArgs, ...)
}

#' Construct a time series of power forecasts, using input training data. Assumes training data already captures
#' differences in magnitude (i.e., power rating) amongst sites. Forecast is NA for times when sun is down.
#'
#' @param x A matrix of training data of dimensions [time x ntrain] for a 1-dimensional forecast
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale Must be "site"
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine' (irrelevant if scale == 'site')
#' @param sun_up_threshold An absolute [MW] threshold on the ensemble members to remove dubious sunrise/sunset valud
#' @param MoreTSArgs An optional dictionary of time-series arguments to the forecast calculation
#' @param ... optional arguments to the prob_forecast object
ts_forecast.matrix <- function(x, start_time, time_step, scale, location, method, sun_up_threshold=0.5, MoreTSArgs=NA, ...) {
  # Check inputs
  if (tolower(scale) !="site"){
    stop("Data and scale mis-match. x must be a 2-dimensional array for 1-dimensional forecasts.")
  }

  new_ts_forecast(x, start_time, time_step, scale, location, method, sun_up_threshold, MoreTSArgs, ...)
}


#' An alternative ts_forecast constructor, which accepts an already calculated list of forecasts.
#' For future work, this could be integrated with copula post-processing with site-level forecasts.
#'
#' @param x A list of prob_forecast objects
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'site', 'region', 'total'
#' @param location A string
#' @param method A string
ts_forecast.list <- function(x, start_time, time_step, scale, location, method) {
  # Check inputs
  if (!any(sapply(x, is.prob_forecast))) stop("Bad input. x must be a list of forecasts.")

  sun_up <- sapply(x, is.prob_forecast)

  dat <- list(start_time = start_time,
              scale = scale,
              location = location,
              time_step = time_step,
              forecasts = x,
              sun_up = sun_up
  )
  structure(dat, class = "ts_forecast")
}

# Helper constructor to calculate forecast from matrix or array data inputs
new_ts_forecast <- function(x, start_time, time_step, scale, location, method,
                                            sun_up_threshold, MoreTSArgs, ...) {
  if (any(names(MoreTSArgs) %in% c("location", "time"))) stop("MoreTSArgs may not include location or time, which are reserved names.")

  sun_up <- apply(x, MARGIN=1, FUN=check_sunup, sun_up_threshold=sun_up_threshold)
  forecasts <- calc_forecasts(x, sun_up, start_time=start_time, time_step=time_step, scale=scale, location=location,
                              method=method, MoreTSArgs=MoreTSArgs, ...)

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
check_sunup <- function(x, sun_up_threshold=0.5){
  return(any(x>=sun_up_threshold, na.rm=TRUE))
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
#' @param MoreTSArgs An optional dictionary of time-series arguments to the forecast calculation
#' @param ... optional arguments to the prob_forecast, including marginal estimator arguments
#' @return A list of forecasts. Forecast is NA for times when sun is down.
calc_forecasts <- function(x, sun_up, start_time, time_step, scale, location, method, MoreTSArgs=NA, ...) {
  class_type <- get_forecast_class(scale, method)
  sub_func <- function(i, time, sun_up){
    if (sun_up) {
      forecast <- tryCatch({
        # If additional time-series arguments are provided, this time points arguments are added to the list of optional arguments
        if (class_type$dim == '1') { # Can't use ifelse function or it will just return the first value
          args <- list(data.input=x[i,], location=location, time=time)
        } else {
          args <- list(data.input=x[i,,], location=location, time=time)
        }

        # Doing this manually to ensure list elements are coerced into a different type
        for (name in names(list(...))) {args[[name]] <- list(...)[[name]]}
        moreargs <- sapply(names(MoreTSArgs), FUN=function(name) return(name=MoreTSArgs[[name]][[i]]), simplify=F)
        for (name in names(moreargs)) {args[[name]] <- moreargs[[name]]}

        return(do.call(class_type$class, args))
      }, error = function(e) {
        e$message <- paste(e$message, "in forecasting time", time, "for location", location)
        stop(e)
      })
      return(forecast)
    } else {return(NA)}}
  return(mapply(sub_func, seq_len(dim(x)[1]), start_time + (seq_len(dim(x)[1])-1)*lubridate::dhours(time_step), sun_up, SIMPLIFY=FALSE))
}

#' Look-up function of the forecast class type
#'
#' @param scale One of 'site', 'region', 'total'
#' @param method One of 'binned' if scale is 'site', else one of gaussian', 'empirical', 'vine'
#' @return A function to initialize a forecast of the desired type
get_forecast_class <- function(scale, method){
  if (tolower(scale) %in% c("site")){
    cls <- switch(tolower(method),
                  "binned" = fc_binned,
                  "kde" = fc_kde,
                  "bma" = fc_bma,
                  "empirical" = fc_empirical,
                  stop(paste('Forecast type', method, 'not recognized for single-site forecasts.', sep=' ')))
    d <- '1'
  } else if (tolower(scale) %in% c("region", "total")){
    cls <- switch(tolower(method),
           "vine" = stop("Not implemented. ts_forecast doesn't have direct handling for vine copula's new list of inputs. "),
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
#' @param window (optional) A vector of (start index, end index) to plot certain time window
plot.ts_forecast <- function(x, tel=NA, window=NA) {
  indices <- get_index_window(x, window)

  probs <- x$forecasts[[min(which(sapply(x$forecasts, FUN=is.prob_forecast)))]]$quantiles$q
  plotdata <- matrix(ncol=length(indices), nrow=length(probs))
  for (i in seq_along(indices)) {
    if (is.prob_forecast(x$forecasts[[indices[i]]])) {
      plotdata[,i] <- x$forecasts[[indices[i]]]$quantiles$x
    } else plotdata[,i] <- 0
  }
  graphics::plot(NULL, xlim=c(0, length(indices)*x$time_step), ylim=c(0, max(plotdata)), xlab="Time [Hrs]", ylab="Power [MW]")
  fanplot::fan(plotdata, data.type='values', probs=probs, fan.col=colorspace::sequential_hcl,
               rlab=NULL, start=x$time_step, frequency=1/x$time_step)
  if (any(!is.na(tel))) {
    actuals_time_step <- x$time_step*length(x)/length(tel)
    graphics::lines(seq(from=actuals_time_step, length.out=length(indices)/actuals_time_step, by=actuals_time_step),
                    tel[((indices[1]-1)*actuals_time_step+1):(indices[length(indices)]*actuals_time_step)], col='chocolate3', lwd=2)
  }
}

#' Export a csv file of [time x quantiles]
#' Note: this follows the convention of e.g., write.csv and write.table, but write isn't an exported generic function
#'
#' @param x A ts_forecast object
#' @param file A file name
#' @param window (optional) A vector of (start index, end index) to plot certain time window
#' @param percentiles (optional -- defaults 1st to 99th percentiles) Select the quantiles associated with which percentiles will be exported
write.ts_forecast <- function(x, file, window=NA, percentiles=1:99) {
  indices <- get_index_window(x, window)

  data <- sapply(percentiles, get_quantile_time_series, x=x, simplify="array")

  colnames(data) <- percentiles
  write.csv(data[indices, ], file=file)
}

# Helper function to select only a subset of the forecast length
get_index_window <- function(x, window) {
  if (all(!(is.na(window))) & length(window)!=2) stop("To use time window, must give a vector of starting and ending indices")

  start <- ifelse(all(is.na(window)), 1, window[1])
  end <- ifelse(all(is.na(window)), length(x), window[2])
  indices <- start:end
}

#' Get a time-series of the forecast value at a particular quantile
#' @param x A ts_forecast object
#' @param alpha Quantile of interest, numeric [0, 100]
#' @return Numeric vector
get_quantile_time_series <- function(x, alpha) {
  eps <- 1e-6
  q <- alpha/100
  timeseries <- sapply(x$forecasts, FUN=function(forecast, q) {
    if (is.prob_forecast(forecast)){v <- forecast$quantiles$x[which(abs(forecast$quantiles$q-q)<eps)]
      if (length(v)==0) stop(paste(alpha,'quantile not available')) else return(v)}
    else return(NA)},
    q=q)
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
#' @param tel A vector of the telemetry values OR a single value of the length of the telemetry values
#' @param fc A vector of data from the time series forecast OR a single value of the length of the forecast values
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @param align Can be "end-of-hour", "half-hour " NaN first, telemetry lags
#' Defaults to "half hour backend". Currently half hour approaches expand forecast the same way.
#' @return list of the telemetry and forecast data vectors, of equal length
equalize_telemetry_forecast_length <- function(tel, fc, agg=TRUE, align="end-of-hour") {
  if (!(align %in% c("end-of-hour", "half-hour"))) stop(paste("Unknown method for align. Given ", align, sep=''))
  if (length(tel) > 1 & length(fc)> 1){
    if ((length(tel) != length(fc)) & (length(tel) < length(fc) | length(tel) %% (2*length(fc)) > 0)) stop("Telemetry length must be equal to or a even multiple of forecast length.")
    fc_length <- length(fc)
    tel_length <- length(tel)
  } else if (length(tel) > 1 & length(fc) == 1) {
    if ((length(tel) != fc) & (length(tel) < fc | length(tel) %% (2*fc) > 0)) stop("Telemetry length must be equal to or a even multiple of forecast length.")
    fc_length <- fc
    tel_length <- length(tel)
  } else if (length(tel) == 1 & length(fc) > 1) {
    if ((tel != length(fc)) & (tel < length(fc) | tel %% (2*length(fc)) > 0)) stop("Telemetry length must be equal to or a even multiple of forecast length.")
    tel_length <- tel
    fc_length <- length(fc)
  } else { # Both tel and fc are singletons
    tel_length <- 1
    fc_length <- 1
  }

  tel_2_fc <- tel_length/fc_length

  # Default translation from telemetry index to forecast index
  index_translation <- function(i) {i}

  if (tel_2_fc > 1) {
    if (align=="end-of-hour") {
      if (agg) {
        tel <- sapply(seq_len(fc_length), function(i) {return(sum(tel[(tel_2_fc*(i-1)+1):(tel_2_fc*i)])/tel_2_fc)})
      } else {
        fc <- rep(fc, each=tel_2_fc)
        # Modified translation function to align more frequent telemetry with forecast:
        index_translation <- purrr::partial(function(i, tel_2_fc) {floor((i-1)/tel_2_fc)+1}, tel_2_fc=tel_2_fc)
      }
    } else { # align == "half-hour"
      if (agg) {
        tel <- c(NaN, sapply(seq_len(fc_length-1), function(i) {return(sum(tel[(tel_2_fc*(i-0.5)+1):(tel_2_fc*(i+0.5))])/tel_2_fc)}))
      } else {
        fc <- c(rep(fc[1], each=tel_2_fc/2), rep(fc[-1], each=tel_2_fc), rep(NA, each=tel_2_fc/2))
        # Modified translation function to align more frequent telemetry with forecast
        index_translation <- purrr::partial(function(i, tel_2_fc) {ceiling((i-0.5*tel_2_fc)/tel_2_fc)+1}, tel_2_fc=tel_2_fc)
      }
    }
  }

  return(list(tel=tel, fc=fc, tel_2_fc=tel_2_fc, translate_forecast_index=index_translation ))
}

#' Get statistics on the prevalence of NaN's in telemetry and sun-up/sun-down discrepancies
#'
#' @param tel A vector of the telemetry values
#' @param len_ts Number of timesteps in the time series forecast
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @param return A list of summary statistics
get_sundown_and_NaN_stats <- function(ts, tel, ...) {
  res <- equalize_telemetry_forecast_length(tel, ts$sun_up, ...)
  sun_up <- res$fc
  tel <- res$tel
  forecast_available <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)$fc

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
              "Telemetry is NaN when sun forecasted down"=tel_sundown_NaN,
              "Sunup missing telemetry rate"=tel_sunup_NaN/(tel_sunup+tel_sunup_NaN),
              "Validatable forecasts"= sum(sun_up & !is.nan(tel) & forecast_available)
         ))
}

#' Get the average estimated CRPS (continuous ranked probability score) for the forecast;
#' the score at each time point is estimated through numerical integration.
#' CRPS characterizes calibration and sharpness together
#'
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param normalize.by (optional) A normalization factor, either a scalar or vector
#' @param ... optional arguments to equalize_telemetry_forecast_length
CRPS_avg <-function(ts, tel, normalize.by=1, ...){
  crps_list <- get_metric_time_series(CRPS, ts, tel, normalize.by, ...)
  return(list(mean=mean(crps_list, na.rm=T), min=min(crps_list, na.rm=T), max=max(crps_list, na.rm=T), sd=stats::sd(crps_list, na.rm=T),
              median=median(crps_list, na.rm=T)))
}

#' Get the quantile-weighted CRPS (continuous ranked probability score) for the forecast
#' Estimated by trapezoidal integration of its quantile decomposition
#' Quantile weighting functions as suggested in Gneiting & Ranjan 2012
#'
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @param quantiles (optional) Sequence of quantiles to integrate over
#' @param qs (optional) A list of quantile scores corresponding to the quantiles
qwCRPS <-function(ts, tel, weighting="none", quantiles=seq(0.01, 0.99, by=0.01), qs=NA){
  if (all(is.na(qs))) {
    if (length(quantiles) != length(qs)) stop("quantile and quantile score vectors must be same length.")
    qs <- QS(ts, tel, quantiles)
  }
  wqs <- weight_QS(qs, quantiles, weighting)
  return(pracma::trapz(quantiles, wqs))
}

#' Get the reliability index (RI) as per Delle Monache 2006
#'
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param nbins Number of bins (original ensemble members) to evaluate
#' @param ... optional arguments to equalize_telemetry_forecast_length
RI <-function(ts, tel, nbins=100, ...){
  counts <- calc_PIT_histogram(ts, tel, nbins, ...)$bin_hits
  return(sum(abs(counts/sum(counts) - 1/nbins)))
}

#' Get the percentile reliability index (PRI) of quantiles 1 to 99.
#' An extension of the Delle Monache 2006 to look at distance from probabilistic calibration (P-P plot).
#' * Note that this looks at the performance at the quantiles, not for each bin -- that is, it is assessed over 99 rather than 100 values
#'
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param ... optional arguments to equalize_telemetry_forecast_length
PRI <-function(ts, tel, ...){
  counts <- calc_PIT_histogram(ts, tel, 100, ...)$bin_hits
  obs_rate <- cumsum(counts/sum(counts))

  return(sum(abs(seq(0.01, 1, by=0.01) - obs_rate)))
}

#' Get a time-series of the forecast performance in terms of one metric
#' @param metric A function for prob_forecast objects
#' @param ts A ts_forecast object
#' @param tel Vector of telemetry
#' @param normalize.by A vector or scalar
#' @param metricArgs (optional) a list of additional arguments to the metrics function
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return A vector
get_metric_time_series <- function(metric, ts, tel, normalize.by, metricArgs=NULL,  ...) {
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc
  normalize.by <- equalize_normalization_factor_length(normalize.by, x$tel, ...)

  metric_list <- sapply(seq_along(forecast_available), function(i) {
    if (forecast_available[i] & !is.na(x$tel[i])) {
      args <- c(list(x=ts$forecasts[[x$translate_forecast_index(i)]], tel=x$tel[i]), metricArgs)
      return(do.call(metric, args)/normalize.by[i])
    } else return(NA)})
}

#' Get weighted quantile score at certain quantile(s)
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param quantiles, numeric [0,1].  Can be a scalar or a vector
#' @return the Quantile score
QS <- function(ts, tel, quantiles) {
  if (length(ts) != length(tel)) stop("Forecast and telemetry must be at same time resolution")
  if (any(quantiles < 0) | any(quantiles > 1)) stop("Quantiles must all be [0,1]")

  thresholds <- t(sapply(100*quantiles, FUN=function(q) {get_quantile_time_series(ts, q)}, simplify="array"))

  valid <- !is.na(ts$forecasts) & !is.na(tel)
  indicator <- sapply(seq_along(quantiles), FUN=function(i) {as.integer(tel[valid] <= thresholds[i, valid])}, simplify="array")

  qs <- sapply(seq_along(quantiles), FUN=function(q) mean(2*(indicator[,q]-quantiles[q])*(thresholds[q,valid]-tel[valid]), na.rm = TRUE))
  return(qs)
}

#' Calculate a vector of weighted quantile scores, emphasizing one or both tails or center
#'
#' @param qs A vector of quantile scores
#' @param quantiles A vector of the quantiles (percentiles) in [0,1]
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @return A vector of the weighted scores
weight_QS <- function(qs, quantiles, weighting="none") {
  if (length(qs) != length(quantiles)) stop("Quantiles and quantile score must be the same length")

  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*quantiles-1)^2,
                    "right" = quantiles^2,
                    "left" = (1-quantiles)^2,
                    "center" = quantiles*(1-quantiles),
                    stop("Weighting method not recognized"))
  return(weights*qs)
}

#' Get Brier score at a certain probability of exceedance
#' This is calculated with a constant probability threshold, rather than a constant value threshold
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param PoE Threshold probability of exceedance, numeric [0,1]
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return the Brier score
Brier_quantile <- function(ts, tel, PoE, ...) {
  if (PoE < 0 | PoE > 1) stop(paste("Probability of exceedance must be [0,1], given ", PoE, '.', sep=''))
  thresholds <- get_quantile_time_series(ts, 100*(1-PoE))

  thresholds <- equalize_telemetry_forecast_length(tel, thresholds, ...)$fc
  res <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- res$fc
  tel <- res$tel

  # Find incidence of telemetry exceeding the threshold of exceedance
  indicator <- as.integer(tel[forecast_available & !is.na(tel)] >= thresholds[forecast_available & !is.na(tel)])
  # Compare Probability of exceedance to its actual incidence
  return(mean((PoE-indicator)^2, na.rm = TRUE))
}

#' Get Brier score at power thresholds
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param thresholds Thresholds in units of the forecast. Can be a scalar or a vector
#' @return the Brier score
Brier <- function(ts, tel, thresholds) {
  if (length(ts) != length(tel)) stop("Forecast and telemetry must be at same time resolution")
  pit <- array(sapply(ts$forecasts, FUN=function(forecast, thresholds) {
    if (is.prob_forecast(forecast)) return(stats::approx(forecast$quantiles$x, forecast$quantiles$q, thresholds, yleft=0, yright=1)$y)
    else return(rep(NA, times=length(thresholds)) )},
    thresholds=thresholds, simplify="array"), dim=c(length(thresholds), length(tel)))

  valid <- !is.na(ts$forecasts) & !is.na(tel)
  indicator <- sapply(seq_along(thresholds), FUN=function(i) as.integer(tel[valid] <= thresholds[i]), simplify="array")

  bs <- sapply(seq_along(thresholds), FUN=function(i) mean((pit[i, valid]-indicator[,i])^2, na.rm = TRUE))

  return(bs)
}

plot_quantile_score <- function(ts, tel, quantiles=seq(0.01, 0.99, by=0.01), weighting="none") {
  qs <- QS(ts, tel, quantiles)
  wqs <- weight_QS(qs, quantiles, weighting)

  if (weighting=="none") label<- "" else label <- paste(weighting,"-weighted", sep="")

  g <- ggplot2::ggplot(data=data.frame(x=quantiles, y=wqs), mapping=ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Quantile") +
    ggplot2::ylab(paste(label, "Quantile Score"))

  plot(g)
}

# Plot Brier score along the quantiles from 1 to 99
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param nmem Number of ensemble members, to illustrate the n+1 bins
plot_brier_by_quantile <- function(ts, tel, nmem = NA) {
  q <- 1:99
  b <- sapply(q, FUN = function(qi) Brier_quantile(ts, tel, (100-qi)/100))

  g <- ggplot2::ggplot(data=data.frame(x=q, y=b), mapping=ggplot2::aes(x=x, y=y)) +
      ggplot2::xlab("Percentile") +
      ggplot2::ylab("Brier score") +
      ggplot2::theme_light()

  if (!is.na(nmem)) {
    xstart <-  seq(0, 100, by=100/(nmem+1)*2)
    xend <- seq(100/(nmem+1), 100, by=100/(nmem+1)*2)
    rects <- data.frame(xstart=xstart[1:length(xend)], xend = xend)
    g <- g + ggplot2::geom_rect(data = rects, mapping=ggplot2::aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), fill = "gray", inherit.aes = FALSE, alpha=0.4)
  }
  g <- g + ggplot2::geom_point() +
        ggplot2::geom_line(data=data.frame(x=q, y=q/100*(100-q)/100))

  plot(g)
}

# Plot Brier score in terms of power
#' @param ts A ts_forecast object
#' @param tel A vector of the telemetry values
#' @param xseq A vector of the power thresholds to use
plot_brier_by_power <- function(ts, tel, xseq) {

  b <- Brier(ts, tel, xseq)

  g <- ggplot2::ggplot(data=data.frame(x=xseq, y=b), mapping=ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_line() +
    ggplot2::xlab("Power [MW]") +
    ggplot2::ylab("Brier score") +
    ggplot2::ggtitle(paste("CRPS:", signif(pracma::trapz(xseq, b), digits = 6), "[MW]"))

  plot(g)
}

equalize_normalization_factor_length <- function(normalize.by, tel, ...) {
  if (length(normalize.by)>1) {
    return(equalize_telemetry_forecast_length(tel, normalize.by, ...)$fc)
  } else {
    return(rep(normalize.by, times=length(tel)))
  }
}

#' Get mean absolute error between the forecast median and the actual value
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param normalize.by (optional) A normalization factor, either a scalar or vector
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return the MAE value
MAE <-function(ts, tel, normalize.by=1, ...) {
  medians <- get_quantile_time_series(ts, 50)
  forecast_available <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)$fc
  x <- equalize_telemetry_forecast_length(tel, medians, ...)
  medians <- x$fc
  normalize.by <- equalize_normalization_factor_length(normalize.by, x$tel, ...)
  indices <- forecast_available & !is.na(x$tel)
  return(mean(abs(medians[indices]-x$tel[indices])/normalize.by[indices], na.rm = TRUE))
}

#' Get average interval score, for an interval from alpha/2 to 1-alpha/2. Negatively oriented (smaller is better)
#' Characterizes sharpness, with a penalty for reliability
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
#' @param normalize.by (optional) A normalization factor, either a scalar or vector
#' @param ... additional optional arguments to equalize_telemetry_forecast_length
#' @return the average IS value
IS_avg <-function(ts, tel, alpha, normalize.by=1, ...) {
  is_list <- get_metric_time_series(IS, ts, tel, normalize.by, metricArgs = list(alpha=alpha), ...)
  return(list(mean=mean(is_list, na.rm=T), min=min(is_list, na.rm=T), max=max(is_list, na.rm=T), sd=stats::sd(is_list, na.rm=T)))
}

#' Get average sharpness, for an interval from alpha/2 to 1-alpha/2. Negatively oriented (smaller is better)
#' Telemetry isn't used in sharpness, but it is only averaged over times when telemetry is available for direct compatibility with other metrics.
#'
#' @param ts A ts_forecast object
#' @param tel A vector of telemetry values.
#' @param alpha Numeric, to identify the (1-alpha)*100% quantile of interest
#' @param normalize.by (optional) A normalization factor, either a scalar or vector
#' @param ... additional optional arguments to equalize_telemetry_forecast_length
#' @return the average sharpness value
sharpness_avg <-function(ts, tel, alpha, normalize.by=1, ...) {
  sharpness_list <- get_metric_time_series(sharpness, ts, tel, normalize.by, metricArgs = list(alpha=alpha), ...)
  return(list(mean=mean(sharpness_list, na.rm=T), min=min(sharpness_list, na.rm=T), max=max(sharpness_list, na.rm=T), sd=stats::sd(sharpness_list, na.rm=T)))
}

#' Plot diagonal line diagram of quantiles + observations (i.e., a P-P plot)
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param ... optional arguments to equalize_telemetry_forecast_length
plot_reliability <- function(ts, tel, ...) {
  x <- get_quantile_reliability(ts, tel, ...)

  ggplot2::ggplot(data.frame(x=x$quantiles, y=x$quantiles), ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data=data.frame(x=x$quantiles[1:(length(x$quantiles)-1)], y=cumsum(x$hit_rate)[1:(length(x$quantiles)-1)]), shape=1) +
    ggplot2::xlab("Nominal Proportion") +
    ggplot2::ylab("Observed Proportion")
}

#' Plot probability integral transform histogram
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param nbins To specify number of bars in the histogram
#' @param ... optional arguments to equalize_telemetry_forecast_length
plot_PIT_histogram <- function(ts, tel, nbins=10, ...) {
  x <- calc_PIT_histogram(ts, tel, nbins, ...)

  ggplot2::ggplot(data.frame(x=x$bin_means, y=x$bin_hits/sum(x$bin_hits)), ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_col(width=x$bin_width) +
    ggplot2::geom_line(data.frame(x=c(0, 1), y=c(1/nbins, 1/nbins)), mapping=ggplot2::aes(x=x, y=y), linetype="dashed") +
    ggplot2::xlab("Probability Integral Transform") +
    ggplot2::ylab("Relative Frequency")
}

#' Calculate data for probability integral transform histogram
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param nbins To specify number of bars in the histogram
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return List of the bin locations, their relative frequency, and the bin width
calc_PIT_histogram <- function(ts, tel, nbins, ...) {
  breaks <- seq(0, 1, length.out =(nbins+1))
  bin_width <- diff(breaks)[1]
  bin_means <- breaks[-1] - bin_width/2
  qr <- get_quantile_reliability(ts, tel, ...)

  # Bin as desired
  bin_rate <- sapply(seq_len(length(breaks)-1), FUN = function(i) sum(qr$hits[qr$quantiles>breaks[i] & qr$quantiles<=breaks[i+1]]))
  return(list(bin_means=bin_means, bin_hits=bin_rate, bin_width=bin_width))
}

#' Evaluate forecast reliability by evaluating the actual hit rate of the quantile bins
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return list of the quantiles and their hit rates
get_quantile_reliability <- function(ts, tel, ...) {
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc

  # Get the list of quantiles that have been evaluated, and add top limit at 1
  quants <- c(ts$forecasts[[min(which(sapply(ts$forecasts, FUN=is.prob_forecast)))]]$quantiles$q, 1)

  # Find time-points where telemetry and forecast data is available
  indices <- which(forecast_available & !is.na(x$tel))
  q_idx_subfunc <- function(i) {
    list_idx <- which(ts$forecasts[[x$translate_forecast_index(i)]]$quantiles$x >= x$tel[i]) # List of quantile indices above telemetry value
    if (length(list_idx) > 0) {idx <- min(list_idx)} else{idx <- length(quants)} # Pick lowest quantile, or 100th percentile if it falls outside distribution
    return(idx)
  }
  q_indices <- sapply(indices, q_idx_subfunc) # Get the indices of the corresponding quantile for each time points
  counts <- sapply(seq_along(quants), function(i) return(sum(q_indices == i)), USE.NAMES = FALSE)

  hit_rate <- counts/length(indices)
  return(list(quantiles=quants, hit_rate=hit_rate, hits=counts))
}
