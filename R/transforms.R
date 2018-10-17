
#' Register generic inverse transform function
#' @param x An object of one of the cdf transform classes
inverse_transform <- function(cdf, x, ...) UseMethod('inverse_transform')


#' Look-up function of transform functions to/from uniform domain
#'
#' @param transform_method 'One of 'empirical', 'kde'
#' @param x Input data
#' @param save (optional) Boolean to skip saving input data for empirical transforms
#' @return list of 'to_uniform' and 'from_uniform' functions
get_transform <- function(transform_method, x, save=TRUE){
  if (!(tolower(transform_method) %in% c("kde", "empirical"))) stop(paste("marginal.transform must be 'kde' or 'empirical'. Given: ", transform_method, sep=''))
  return(switch(tolower(transform_method),
                          "kde" = cdf_kde(x),
                          "empirical" = if (save) {cdf_em(x)} else {cdf_em(NULL)}))
}

# -------------------------------------------------------------------------------------------------------------

#' A structure for transforming back and forth to uniform domain using a kernel density estimate (KDE) of the input data, truncated at 0
#' Using KernSmooth based on https://vita.had.co.nz/papers/density-estimation.pdf
#'
#' @param x A vector of samples
#' @param bandwidth Kernel bandwidth (optional -- otherwise uses oversmoothed bandwidth selector)
#' @param kernel Smoothing kernel for KernSmooth bkde
#' @return
cdf_kde <- function(x, ...,
                    bandwidth=NA, kernel='norm') {

  if (is.na(bandwidth)) {
    pdf <- KernSmooth::bkde(x, kernel=kernel, range.x=c(0, 1.1*max(x)),  truncate=TRUE)
  } else
    pdf <- KernSmooth::bkde(x, kernel=kernel, bandwidth=bandwidth, range.x=c(0, 1.1*max(x)),  truncate=TRUE)
  # Calculate cumulative distribution function with zero order hold
  cdf_raw <- cumsum(pdf$y)*(pdf$x[2]-pdf$x[1])
  # Correct to ensure CDF maximum is 1
  cdf <- cdf_raw/max(cdf_raw)
  dat <- list('x'=pdf$x,
              'pdf'=pdf$y,
              'cdf'=cdf)
  x <- structure(dat, class = c("cdf_kde"))
  return(x)
}

#' Check probabilistic forecast class
is.cdf_kde <- function(x) inherits(x, "cdf_kde")

#' Inverse transform from uniform to variable domain
#'
#' @param cdf_kde A cdf_kde object
#' @param u A vector of evaluation points to transform
#' @return A vector linearly interpolated from the uniform to variable domain
inverse_transform.cdf_kde <- function(cdf_kde, u) {
 return(approx(cdf_kde$cdf, cdf_kde$x, xout=u, rule=2)$y)
}

#' Transform from variable to uniform domain
#'
#' @param cdf A cdf_kde object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
transform.cdf_kde <- function(cdf, x) {
  return(approx(cdf$x, cdf$cdf, xout=x, rule=2)$y)
}

plot.cdf_kde <- function(cdf) {
  plot(cdf$x, cdf$cdf, main="KDE CDF", xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")),
       ylab='Uniform domain U=F(X)')
}

plot_pdf <- function(cdf) {
  plot(cdf$x, cdf$pdf, main="KDE PDF", xlab="Variable domain X", ylab='Relative frequency')
}

min_x <- function(cdf){
  return(min(cdf$x))
}

max_x <- function(cdf){
  return(max(cdf$x))
}

# -----------------------------------------------------------------------------------------------------

#' A structure for transforming back and forth to uniform domain using an empirical CDF of the input data
#'
#' @param x A vector of samples
#' @return
cdf_em <- function(x, ...) {
  x <- structure(list('samples'=x), class = c("cdf_em"))
  return(x)
}

#' Check probabilistic forecast class
is.cdf_em <- function(x) inherits(x, "cdf_em")

#' Inverse transform from uniform to variable domain
#'
#' @param cdf_em A cdf_em object
#' @param u A vector of evaluation points to transform
#' @return A vector linearly interpolated from the uniform to variable domain
inverse_transform.cdf_em <- function(cdf_em, u) {
  return(stats::quantile(cdf_em$samples, u, type=4, names=FALSE))
}

#' Transform from variable to uniform domain
#'
#' @param cdf_em A cdf_em object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
transform.cdf_em <- function(cdf_em, x) {
  return(rvinecopulib::pseudo_obs(x))
}

plot.cdf_em <- function(cdf_em) {
  plot(ecdf(cdf_em$samples), main="Empirical CDF", xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")),
       ylab='Uniform domain U=F(X)')
}
