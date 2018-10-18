#'
#' Assumes that the forecast data and historical analogs are at the same time-resolution
#'
#' @param f_test matrix of the forecast data to fit on, rows for time-steps, columns for physical variables
#' @param h_train matrix of historical forecast data, rows for time-steps, columns for physical variables
#' @param h_real Historical realized value of interest, e.g., power (equivalent to a kNN classification)
#' @param n Integer, number of historical analogs to pick
#' @param features Vector of the column names to use for similarity measurement
#' @param weights Vector of weights to use for each feature
#' @return A matrix of the analogs, one analog per row, columns for physical variables + realized value
get_historical_analogs <- function(f_test, h_train, h_real, n, features, weights) {
  # Error check arguments
  if (dim(f_test)[1] %% 2 == 0) stop('Forecast data must have an odd number of time-steps for centerd analog matching.')
  if (length(features) != length(weights)) stop(paste('Given', length(features), 'features and', length(weights), 'weights. Must be same length.', sep=' '))
  if (sum(weights) != 1) stop('Weights must sum to 1.')
  if (n < 1) stop(paste('At least 1 analog required. Requested', n, 'analogs.', sep=' '))

  # How does this work with h_train netCDF format?
  sigmas <- apply(h_train[,features], 2, sd)

  metrics <- vapply(1:dim(h_train)[1], delle_monache_distance, numeric(1),
                     f=f_test, h=h_train, features=features, weights=weights, sigmas=sigmas)

  analog_metrics <- metrics[order(metrics, na.last = TRUE)[1:n]]
  analogs <- h_real[order(metrics)[1:n]]
  forecasts <- h_train[order(metrics)[1:n],]

  return(list('real'=analogs, 'distance'=analog_metrics, 'forecast'=forecasts))
}

#' Calculate Delle Monache distance metric for the physical features for a potential analog
#'
#' @param t_prime Centered time index of the potential analog
#' @param f forecast of features the over the time interval of interest
#' @param h historical forecasts of the features
#' @param features Vector of the column names to use for similarity measurement
#' @param weights Vector of weights to use for each feature
#' @param sigmas Standard deviations of the features over the historical record
#' @return Delle Monache distance
delle_monache_distance <- function(t_prime, f, h, features, weights, sigmas) {
  half_interval <- (dim(f)[1] - 1)/2

  if (t_prime <= half_interval | t_prime >= dim(h)[1]-half_interval+1) {
    return(NA)
  } else {
   feat_dist <- mapply(feature_distance, features, weights, sigmas, MoreArgs = list(f=f, h=h[(t_prime-half_interval):(t_prime+half_interval), features]))
   return(sum(feat_dist))
  }
}

#' Calculate the distance metric for one physical feature for a potential analog
#'
#' @param feature a column name
#' @param weight weight [0,1] assigned to that feature
#' @param sigma Standard deviation of the feature over the historical record
#' @param f forecast of features the over the time interval of interest
#' @param h possible historical analog of the features, sliced to equivalent time interval
#' @return feature distance
feature_distance <- function(feature, weight, sigma, f, h) {
 return(weight/sigma*sqrt(sum((f[,feature]-h[,feature])^2)))
}
