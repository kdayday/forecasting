#' A structure for transforming back and forth to uniform domain
#'
#' @param x A vector of samples
#' @param method Method for estimating the marginal distribution. One of 'empirical', 'kde1d', 'geenens', 'kernsmooth'
#' @param ... Optional arguments to the distribution estimator. 'xmax' can be used for 'kde1d' and 'geenens'
#' @return
marg_transform <- function(x, method='geenens', ... ) {
  if (class(method) != 'character') stop('Method selection must be a character name.')

  func <- kde_lookup(method)

  # Get selected marginal estimate
  res <- func(x, ...)

  if (!all(is.nan(res$d)) & any(res$d < 0)) warning("Negative probability density values.") # Ignores empirical estimate

  dat <- list('x'=res$x,
              'd'=res$d,
              'u'=res$u,
              'xmin' = min(res$x),
              'xmax' = max(res$x),
              'method' = method)
  if (dat$xmin != 0) warning(paste("Lower boundary of support is ", dat$xmin, " rather than 0.", sep=''))

  x <- structure(dat, class = c("marg_transform"))
  return(x)
}

#' Check marginal transform class
is.marg_transform <- function(x) inherits(x, "marg_transform")

summary.marg_transform <- function(x) {
  print("Marginal transform")
  print(paste("Method:", x$method, sep=' '))
  print(paste("Estimated over support from 0 to", x$xmax, sep=' '))
  print(paste("Max CDF value:", max(x$u), sep=' '))
}

#' Plot a marginal transform
#' @param c A marg_transform object
plot.marg_transform <- function(c) {
  graphics::plot(c$x, c$u, xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")), type='l', lwd=2,
       ylab='Uniform domain U=F(X)')
}

#' Plot a marg_transform PDF
#'
#' @param c A marg_transform object
plot_pdf <- function(c, col='black') {
  if (all(is.nan(c$d))) stop('Probability density estimate unavailable. (Perhaps using empirical estimate?)')
  graphics::plot(c$x, c$d, xlab="Variable domain X", ylab='Probability density', type='l', lwd=2, col=col)
}

# ---------------------------------------------------------
# Transformation functions

#' Inverse transform from uniform to variable domain
#'
#' @param c A marg_transform object
#' @param u A vector of evaluation points to transform
#' @return A vector linearly interpolated from the uniform to variable domain
from_uniform <- function(c, u) {
  if (any(u<=0 | u>=1)) stop("Evaluation point(s) for uniform transform must be in (0,1).")
  return(stats::approx(c$u, c$x, xout=u, rule=2)$y)
}

#' Transform from variable to uniform domain
#'
#' @param c A marg_transform object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to uniform domain
to_uniform <- function(c, x) {
  if (any(x<c$xmin | x >c$xmax)) warning("Evaluation point(s) beyond the variable range of the CDF.  ")
  return(stats::approx(c$x, c$u, xout=x, yleft = 0, yright = 1)$y)
}

#' Transform from variable to probability density
#'
#' @param c A marg_transform object
#' @param x A vector of evaluation points to transform
#' @return A vector linearly interpolated from the variable to probability density
to_probability <- function(c, x) {
  if (any(x<c$xmin | x >c$xmax)) warning("Evaluation point(s) beyond the variable range of the PDF.  ")
  if (is.nan(c$d)) stop("Probability vector is NaN. (Undefined for empirical transform).")
  return(stats::approx(c$x, c$d, xout=x, yleft = 0, yright = 0)$y)
}
