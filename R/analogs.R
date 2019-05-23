#'
#' Assumes that the forecast data and historical analogs are at the same time-resolution
#'
#' @param f_test matrix of the forecast data to fit on, rows for time-steps, columns for physical variables
#' @param h_train matrix of historical forecast data, rows for time-steps, columns for physical variables
#' @param h_real Historical realized value of interest, e.g., power (equivalent to a kNN classification)
#' @param n Integer, number of historical analogs to pick
#' @param weights Vector of weights to use for each feature
#' @return A matrix of the analogs, one analog per row, columns for physical variables + realized value
get_historical_analogs <- function(f_test, h_train, h_real, n, weights) {
    # Error check arguments
  if (dim(f_test)[1] %% 2 == 0) stop('Forecast data must have an odd number of time-steps for centered analog matching.')
  if (dim(f_test)[2] != length(weights)) stop("Must have same number of weights as physical features")
  if (sum(weights) != 1) stop('Weights must sum to 1.')
  if (n < 1) stop(paste('At least 1 analog required. Requested', n, 'analogs.', sep=' '))

  # How does this work with h_train netCDF format?
  sigmas <- apply(h_train, 2, sd, na.rm=T)

  metrics <- vapply(1:dim(h_train)[1], delle_monache_distance, numeric(1),
                     f=f_test, h=h_train, weights=weights, sigmas=sigmas)

  # NA's get removed
  indices <- order(metrics, na.last=NA)
  # If there are fewer available data points than requested, use them all; otherwise, grab best n fits
  if (length(indices) > n)
    indices <- indices[1:n]

  analog_metrics <- metrics[indices]
  analogs <- h_real[indices]
  forecasts <- h_train[indices,]

  return(list('real'=analogs, 'distance'=analog_metrics, 'forecast'=forecasts))
}

#' Calculate Delle Monache distance metric for the physical features for a potential analog
#'
#' @param t_prime Centered time index of the potential analog
#' @param f forecast of features the over the time interval of interest
#' @param h historical forecasts of features
#' @param weights Vector of weights to use for each feature
#' @param sigmas Standard deviations of the features over the historical record
#' @return Delle Monache distance
delle_monache_distance <- function(t_prime, f, h, weights, sigmas) {
  half_interval <- (dim(f)[1] - 1)/2

  if (t_prime <= half_interval | t_prime >= dim(h)[1]-half_interval+1) {
    return(NA)
  } else {
   feat_dist <- mapply(feature_distance, seq_len(dim(f)[2]), weights, sigmas, MoreArgs = list(f=f, h=h[(t_prime-half_interval):(t_prime+half_interval),]))
   return(sum(feat_dist))
  }
}

#' Calculate the distance metric for one physical feature for a potential analog
#' Raw distance is 0 if BOTH forecast and historical are missing
#' to accommodate windows that overlap sunrise/set. Distance is still NA if only one is missing.
#'
#' @param i index of feature column
#' @param weight weight [0,1] assigned to this feature
#' @param sigma Standard deviation of the feature over the historical record
#' @param f forecast of features the over the time interval of interest
#' @param h possible historical analog of the features, sliced to equivalent time interval
#' @return feature distance
feature_distance <- function(i, weight, sigma, f, h) {
  diff <- f[,i]-h[,i]
  # Overwrite with distance of 0 if both the forecast and historical were NA -- for times overlapping with sunrise/set
  diff[is.na(f[,i]) & is.na(h[,i])] <- 0
 return(weight/sigma*sqrt(sum((diff)^2)))
}


