#' Initialize a time series of power forecasts. Assumes training data already captures
#' differences in magnitude (i.e., power rating) amongst sites. Forecast is NA for times when sun is down.
#'
#' @param x A matrix or array of training data of dimensions [time x ntrain x nsites] for n-dimensional forecast or [time x ntrain] for a 1-dimensional forecast
#' @param start_time A lubridate time stamp
#' @param time_step Time step in hours
#' @param scale One of 'site', 'region', 'total'
#' @param location A string
#' @param method One of 'gaussian', 'empirical', 'vine' (irrelevant if scale == 'site')
#' @param sun_up_threshold An absolute [MW] threshold on the ensemble members to remove dubious sunrise/sunset valud
#' @param MoreTSArgs An optional dictionary of time-series arguments to the forecast calculation
#' @param ... optional arguments to the prob_forecast object
ts_forecast <- function(x, start_time, time_step, scale, location, method, sun_up_threshold=0.5, MoreTSArgs=NA, ...) {
  # Check inputs
  if (!is.array(x)) stop("Bad input. x must be an array.")
  if (!((length(dim(x))==2 & tolower(scale) %in% c("site", "s")) | (length(dim(x))==3 & tolower(scale) %in% c("region", "total", "r", "t")))){
    stop("Data and scale mis-match. x must be a 2-dimensional array for 1-dimensional forecasts or a 3-dimensional array for n-dimensional forecasts.")
  }
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
#' @param method One of 'rank' if scale is 'site', else one of gaussian', 'empirical', 'vine'
#' @return A function to initialize a forecast of the desired type
get_forecast_class <- function(scale, method){
  if (tolower(scale) %in% c("site", "s")){
    cls <- switch(tolower(method),
                  "rank" = prob_1d_rank_forecast,
                  "r" = prob_1d_rank_forecast,
                  "kde" = prob_1d_kde_forecast,
                  "k" = prob_1d_kde_forecast,
                  "bma" = prob_1d_bma_forecast,
                  "b" = prob_1d_bma_forecast,
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
  probs <- x$forecasts[[min(which(sapply(x$forecasts, FUN=is.prob_forecast)))]]$quantiles$q
  plotdata <- matrix(ncol=length(x), nrow=length(probs))
  for (i in seq_along(x$forecasts)) {
    if (is.prob_forecast(x$forecasts[[i]])) {
      plotdata[,i] <- x$forecasts[[i]]$quantiles$x
    } else plotdata[,i] <- 0
  }
  graphics::plot(NULL, xlim=c(0, length(x)*x$time_step), ylim=c(0, max(plotdata)), xlab="Time [Hrs]", ylab="Power [MW]")
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
#' @param tel A vector of the telemetry values
#' @param fc A vector of data from the time series forecast
#' @param agg Boolean, TRUE to aggregate telemetry to forecast resolution, FALSE to expand forecast to telemetry resolution
#' @param align Can be "top of hour", "half hour backend" (NaN first, telemetry lags), or "half hour frontend" (NaN last, telemetry leads).
#' Defaults to "half hour backend". Currently half hour approaches expand forecast the same way.
#' @return list of the telemetry and forecast data vectors, of equal length
equalize_telemetry_forecast_length <- function(tel, fc, agg=TRUE, align="end-of-hour") {
  if (!(align %in% c("end-of-hour", "half-hour-backend", "half-hour-frontend"))) stop(paste("Unknown method for align. Given ", align, sep=''))
  if (length(tel) != length(fc) & (length(tel) < length(fc) | length(tel) %% (2*length(fc)) > 0)) stop("Telemetry length must be equal to or a even multiple of forecast length.")

  tel_2_fc <- length(tel)/length(fc)

  # Default translation from telemetry index to forecast index
  index_translation <- function(i) {i}

  if (tel_2_fc > 1) {
    if (align=="end-of-hour") {
      if (agg) {
        tel <- sapply(seq_along(fc), function(i) {return(sum(tel[(tel_2_fc*(i-1)+1):(tel_2_fc*i)])/tel_2_fc)})
      } else {
        fc <- rep(fc, each=tel_2_fc)
        # Modified translation function to align more frequent telemetry with forecast:
        index_translation <- purrr::partial(function(i, tel_2_fc) {floor((i-1)/tel_2_fc)+1}, tel_2_fc=tel_2_fc)
      }
    } else { # align == "half hour" of some variety
      if (agg) {
        if (align=="half-hour-backend"){
          tel <- c(NaN, sapply(seq_along(fc[-1]), function(i) {return(sum(tel[(tel_2_fc*(i-0.5)+1):(tel_2_fc*(i+0.5))])/tel_2_fc)}))
        } else {# align == "half hour frontend"
          tel <- c(sapply(seq_along(fc[-1]), function(i) {return(sum(tel[(tel_2_fc*(i-0.5)+1):(tel_2_fc*(i+0.5))])/tel_2_fc)}), NaN)
        }
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
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc
  normalize.by <- equalize_normalization_factor_length(normalize.by, x$tel, ...)
  crps_list <- sapply(which(forecast_available & !is.na(x$tel)), function(i) return(CRPS(ts$forecasts[[x$translate_forecast_index(i)]], x$tel[i])/normalize.by[i]))
  return(list(mean=mean(crps_list), min=min(crps_list), max=max(crps_list), sd=stats::sd(crps_list)))
}

#' Get Brier score at a certain probability of exceedance
#' This is calculated with a constant probability threshold, rather than a constant value threshold
#'
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param PoE Threshold probability of exceedance, numeric [0,1]
#' @param ... optional arguments to equalize_telemetry_forecast_length
#' @return the Brier score
Brier <- function(ts, tel, PoE, ...) {
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

# Plot Brier score along the quantiles from 1 to 99
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param nmem Number of ensemble members, to illustrate the n+1 bins
plot_brier <- function(ts, tel, nmem = NA) {
  q <- 1:99
  b <- sapply(q, FUN = function(qi) Brier(ts, tel, (100-qi)/100))

  g <- ggplot2::ggplot(data=data.frame(x=q, y=b), mapping=ggplot2::aes(x=x, y=y)) +
      ggplot2::xlab("Quantile") +
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
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc
  normalize.by <- equalize_normalization_factor_length(normalize.by, x$tel, ...)
  is_list <- sapply(which(forecast_available & !is.na(x$tel)), function(i) {IS(ts$forecasts[[x$translate_forecast_index(i)]], x$tel[i], alpha=alpha)/normalize.by[i]})
  return(list(mean=mean(is_list), min=min(is_list), max=max(is_list), sd=stats::sd(is_list)))
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
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc
  normalize.by <- equalize_normalization_factor_length(normalize.by, x$tel, ...)
  sharpness_list <- sapply(which(forecast_available & !is.na(x$tel)), function(i) {sharpness(ts$forecasts[[i]], alpha=alpha)/normalize.by[i]})
  return(list(mean=mean(sharpness_list), min=min(sharpness_list), max=max(sharpness_list), sd=stats::sd(sharpness_list)))
}

#' Plot diagonal line diagram of quantiles + observations
#' @param ts A ts_forecast object
#' @param tel A list of the telemetry values
#' @param ... optional arguments to equalize_telemetry_forecast_length
plot_reliability <- function(ts, tel, ...) {
  x <- get_quantile_reliability(ts, tel, ...)

  graphics::plot(x$quantiles, x$quantiles, type="l", lty=2, xlab="Nominal",
                 ylab="Observed")
  graphics::lines(x$quantiles[1:(length(x$quantiles)-1)], cumsum(x$hit_rate)[1:(length(x$quantiles)-1)], type='b', lty=1, pch=1)
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
#' @param quants Provide a vector of bin breakpoints, or NA to use the given quantiles
#' @return list of the quantiles and their hit rates
get_quantile_reliability <- function(ts, tel, ..., quants=NA) {
  x <- equalize_telemetry_forecast_length(tel, !is.na(ts$forecasts), ...)
  forecast_available <- x$fc

  # Get the list of quantiles that have been evaluated, and add top limit at 1
  if (all(is.na(quants))) {
    quants <- c(ts$forecasts[[min(which(sapply(ts$forecasts, FUN=is.prob_forecast)))]]$quantiles$q, 1)
  }

  # Find time-points where telemetry and forecast data is available
  indices <- which(forecast_available & !is.na(x$tel))
  q_idx_subfunc <- function(i) {
    list_idx <- which(ts$forecasts[[x$translate_forecast_index(i)]]$quantiles$x > x$tel[i]) # List of quantile indices above telemetry value
    if (length(list_idx) > 0) {idx <- min(list_idx)} else{idx <- length(quants)} # Pick lowest quantile, or 100th percentile if it falls outside distribution
    return(idx)
  }
  q_indices <- sapply(indices, q_idx_subfunc) # Get the indices of the corresponding quantile for each time points
  counts <- sapply(seq_along(quants), function(i) return(sum(q_indices == i)), USE.NAMES = FALSE)

  hit_rate <- counts/length(indices)
  return(list(quantiles=quants, hit_rate=hit_rate, hits=counts))
}
