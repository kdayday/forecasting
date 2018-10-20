
#' Register generic variable-to-uniform domain transform function
#' @param x An object of one of the cdf transform classes
to_uniform <- function(cdf, x, ...) UseMethod('to_uniform')

#' Register generic uniform-to-variable domain transform function
#' @param x An object of one of the cdf transform classes
from_uniform <- function(cdf, x, ...) UseMethod('from_uniform')

#' Register generic variable-to-probability density transform function
#' @param x An object of one of the cdf transform classes
to_probability <- function(cdf, x, ...) UseMethod('to_probability')

#' Look-up function of transform functions to/from uniform domain
#'
#' @param x Input data
#' @param transform_method 'One of 'empirical', 'kde', default to 'kde'
#' @return list of 'to_uniform' and 'from_uniform' functions
get_transform <- function(x, transform_method='kde'){
  if (!(tolower(transform_method) %in% c("kde", "empirical"))) stop(paste("marginal.transform must be 'kde' or 'empirical'. Given: ", transform_method, sep=''))
  if (tolower(transform_method) == 'empirical') {warning("Empirical transform functionality is incomplete.")}
  return(switch(tolower(transform_method),
                          "kde" = cdf_kde(x),
                          "empirical" = cdf_em(x)))
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
              'cdf'=cdf,
              'xmin' = min(pdf$x),
              'xmax' = max(pdf$x))
  x <- structure(dat, class = c("cdf_kde"))
  return(x)
}

#' Check probabilistic forecast class
is.cdf_kde <- function(x) inherits(x, "cdf_kde")

#' Inverse transform from uniform to variable domain
#'
#' @param c A cdf_kde object
#' @param u A vector of evaluation points to transform
#' @return A vector linearly interpolated from the uniform to variable domain
from_uniform.cdf_kde <- function(c, u) {
  if (any(u<0 | u>1)) stop("Evaluation point(s) for uniform transform must be in [0,1].")
 return(approx(c$cdf, c$x, xout=u, rule=2)$y)
}

#' Transform from variable to uniform domain
#'
#' @param c A cdf_kde object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
to_uniform.cdf_kde <- function(c, x) {
  if (any(x<c$xmin | x >c$xmax)) warning("Evaluation point(s) beyond the variable range of the CDF.  ")
  return(approx(c$x, c$cdf, xout=x, yleft = 0, yright = 1)$y)
}

#' Transform from variable to probability density
#'
#' @param c A cdf_kde object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to probability density
to_probability.cdf_kde <- function(c, x) {
  if (any(x<c$xmin | x >c$xmax)) warning("Evaluation point(s) beyond the variable range of the PDF.  ")
  return(approx(c$x, c$pdf, xout=x, yleft = 0, yright = 0)$y)
}

#' Plot a KDE CDF transform
#'
#' @param c A cdf_kde object
plot.cdf_kde <- function(c) {
  plot(c$x, c$cdf, main="KDE CDF", xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")),
       ylab='Uniform domain U=F(X)')
}

#' Plot a KDE PDFobject
#'
#' @param c A cdf_kde
plot_pdf <- function(cdf) {
  plot(cdf$x, cdf$pdf, main="KDE PDF", xlab="Variable domain X", ylab='Relative frequency')
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
from_uniform.cdf_em <- function(cdf_em, u) {
  return(stats::quantile(cdf_em$samples, u, type=4, names=FALSE))
}

#' Transform from variable to uniform domain
#'
#' @param cdf_em A cdf_em object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
to_uniform.cdf_em <- function(cdf_em, x) {
  return(rvinecopulib::pseudo_obs(x))
}

#' Transform from variable to uniform domain
#'
#' @param cdf_em A cdf_em object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
to_probability.cdf_em <- function(cdf_em, x) {
  stop('Not implemented. Just use KDE of probability density.')
}

#' Plot an empirical CDF transform
#'
#' @param c A cdf_em object
plot.cdf_em <- function(cdf_em) {
  plot(ecdf(cdf_em$samples), main="Empirical CDF", xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")),
       ylab='Uniform domain U=F(X)')
}
