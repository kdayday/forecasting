#'
#' Assumes that the forecast data and historical analogs are at the same time-resolution
#' To calculate standard deviation of each feature, zero values are skipped.
#'
#' @param f_test matrix of the forecast data to fit on [time along matching window x physical feature]
#' @param h_train array of historical forecast data [potential analog time x time along matching window x physical feature]. Matching window time can be a singleton dimension
#' @param h_real Historical realized value of interest, e.g., power (equivalent to a kNN classification)
#' @param n Integer, number of historical analogs to pick
#' @param weights Vector of weights to use for each feature
#' @param sigmas (optional) Vector of standard deviations to use for each feature. Can be re-calculated if feature has redundancies to make matrix structure.
#' @return A list of the analogs, including observed value, the forecast along the matching window, and the distance metric.
get_historical_analogs <- function(f_test, h_train, h_real, n, weights, sigmas=FALSE) {
    # Error check arguments
  matching_time <- dim(f_test)[1]
  n_features <- dim(f_test)[2]
  if (matching_time %% 2 == 0) stop('Forecast data must have an odd number of time-steps for centered analog matching.')
  if (matching_time != dim(h_train)[2]) stop('Forecast data and test data must have matching windows of the same length.')
  if (n_features != length(weights)) stop("Must have same number of weights as physical features")
  if (dim(h_train)[1] != length(h_real)) stop("Historical forecast and telemetry data must be the length, i.e., already at the same time resolution.")
  if (sum(weights) != 1) stop('Weights must sum to 1.')
  if (n < 1) stop(paste('At least 1 analog required. Requested', n, 'analogs.', sep=' '))



  # There should be the same number of standard deviations as weights
  if (!sigmas) sigmas <- apply(h_train, 3, FUN=function(i) sd(i[i>0], na.rm=T))

  metrics <- vapply(seq_along(h_real), function(i, ...) return(ifelse(is.na(h_real[i]), NA, delle_monache_distance(h_train[i,,], ...))),
                    FUN.VALUE=numeric(1), f=f_test, weights=weights, sigmas=sigmas)

  # NA's get removed
  indices <- order(metrics, na.last=NA)
  # If there are no viable analogs, throw error
  if (length(indices) == 0)
    stop("No viable analogs found")

  # If there are fewer available data points than requested, use them all and top up with NaN's for remainder; otherwise, grab best n fits
  if (length(indices) >= n) {
    analog_metrics <- metrics[indices[1:n]]
    analogs <- h_real[indices[1:n]]
    forecasts <- h_train[indices[1:n],,]
  } else {
    n_missing <- n-length(indices)
    analog_metrics <- c(metrics[indices], rep(NaN, times=n_missing))
    analogs <- c(h_real[indices], rep(NaN, times=n_missing))
    forecasts <- aperm(array(c(h_train[indices,,], rep(NaN, times=n_missing*prod(dim(f_test)))), dim=c(matching_time, n_features, n)), c(3, 1, 2))
  }

  return(list('obs'=analogs, 'distance'=analog_metrics, 'forecast'=forecasts))
}

#' Calculate Delle Monache distance metric for the physical features for a potential analog
#'
#' @param h a potential analog, [time along matching window x physical feature]
#' @param f forecast of features, [time along matching window x physical feature]
#' @param weights Vector of weights to use for each feature
#' @param sigmas Standard deviations of the features over the historical record
#' @return Delle Monache distance
delle_monache_distance <- function(h, f, weights, sigmas) {
  feat_dist <- sapply(seq_len(dim(f)[2]), FUN=function(i) return(feature_distance(weights[i], sigmas[i], f[,i], h[,i])))
  return(sum(feat_dist))
}

#' Calculate the distance metric for one physical feature for a potential analog
#' Raw distance is 0 if BOTH forecast and historical are missing
#' to accommodate windows that overlap sunrise/set. Distance is still NA if only one is missing.
#' Matrix slicing is done internally to handle half_interval of 0 or matrix of single feature (to avoid vectors)
#'
#' @param weight weight [0,1] assigned to this feature
#' @param sigma Standard deviation of the feature over the historical record
#' @param f A vector forecast of one features
#' @param h An potential analog forecast of one feature
#' @return feature distance
feature_distance <- function(weight, sigma, f, h) {
  diff <- f-h
  # Overwrite with distance of 0 if both the forecast and historical were NA -- for times overlapping with sunrise/set
  diff[is.na(f) & is.na(h)] <- 0
 return(weight/sigma*sqrt(sum((diff)^2)))
}
