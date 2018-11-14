#' A structure for transforming back and forth to uniform domain
#'
#' @param x A vector of samples
#' @param method Method for estimating the marginal distribution. One of 'empirical', 'kde1d', 'geenens', 'kernsmooth'
#' @param ... Optional arguments to the distribution estimator. 'xmin' and 'xmax' can be used for 'kde1d' and 'geenens'
#' @return
marg_transform <- function(x, method='geenens', ... ) {
  if (class(method) != 'character') stop('Method selection must be a character name.')

  func <- switch(method,
         empirical = probempirical,
         kde1d = probkde1d,
         geenens = probtranskde,
         kernsmooth = probkde,
         stop(paste("Marginal distribution estimator ", method, " not recognized.", sep='')))

  # Get selected marginal estimate
  res <- func(x, ...)

  if (!all(is.na(res$d)) & any(res$d < 0)) warning("Negative probability density values.") # Ignores empirical estimate

  dat <- list('x'=res$x,
              'd'=res$d,
              'u'=res$u,
              'xmin' = min(res$x),
              'xmax' = max(res$x),
              'method' = method)

  x <- structure(dat, class = c("marg_transform"))
  return(x)
}

#' Check marginal transform class
is.marg_transform <- function(x) inherits(x, "marg_transform")

summary.marg_transform <- function(x) {
  print("Marginal transform")
  print(paste("Method:", x$method, sep=' '))
  print(paste("Estimated over support from", xmin, 'to', xmax, sep=' '))
  print(paste("Area under CDF:", pracma::trapz(x$x, x$u, sep=' ')))
}

#' Plot a marginal transform
#' @param c A marg_transform object
plot.marg_transform <- function(c) {
  graphics::plot(c$x, c$u, xlab=expression(paste("Variable domain X=", F^-1, "(", U, ")", sep = "")),
       ylab='Uniform domain U=F(X)')
}

#' Plot a marg_transform PDF
#'
#' @param c A marg_transform object
plot_pdf <- function(c) {
  if (all(is.na(c$d))) stop('Probability density estimate unavailable. (Perhaps using empirical estimate?)')
  graphics::plot(c$x, c$d, xlab="Variable domain X", ylab='Probability density')
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
  if (is.na(c$d)) stop("Probability vector is NA. (Undefined for empirical transform).")
  return(stats::approx(c$x, c$d, xout=x, yleft = 0, yright = 0)$y)
}

# -------------------------------------------------------------------------
# Containers for existing distribution estimation functions

#' Get continuous KDE using KernSmooth based on https://vita.had.co.nz/papers/density-estimation.pdf
#'
#' @param x A vector of samples
#' @param ... Optional inputs to the bkde function
#' @return A list of the evaluation points, density, and cumulative distribution
probkde <- function(x, ...) {
  pdf <- KernSmooth::bkde(x, ...)
  cdf <- pracma::cumtrapz(pdf$x, pdf$y)
  return(list(x=pdf$x, d=as.numeric(pdf$y), u=as.numeric(cdf)))
}

#' Get KDE using kde1d
#'
#' @param x A vector of samples
#' @param ... Optional inputs to the bkde function
#' @return A list of the evaluation points, density, and cumulative distribution
probkde1d <- function(x, ...) {
  stop("Not implemented yet, undetermined bugs.")
  pdf <- kde1d::kde1d(x, ...)
  cdf <- pracma::cumtrapz(pdf$grid_points, pdf$values)
  return(list(x=pdf$grid_points, d=as.numeric(pdf$values), u=as.numeric(cdf)))
}

#' Get probability using empirical CDF, with a default minimum x value of 0
#'
#' @param x A vector of samples
#' @param xmax maximum
#' @return A list of the evaluation points, density, and cumulative distribution
probempirical <- function(x, xmax=max(x)) {

  # Pick the largest maximum value, in case actual data exceeds rating
  xmax = max(xmax, max(x))
  if (xmax > max(x)) {
    xseq <- c(0, sort(x), xmax)
  } else {
      xseq <- c(0, sort(x))
  }

  # Aggregate to combine duplicate values
  agg <- stats::aggregate(xseq, list(value = xseq), length)
  xout <- as.numeric(unlist(agg['value']))
  u <- as.numeric(unlist((cumsum(agg['x'])-1)/(length(xseq)-1)))

  return(list(x=xout, u=u, d=NA))
}

# ----------------------------------------------------------------------------
# Functions for bounded support transformations with nearest neighbors bandwidths

# Author: Gery Geenens
myintegral<-function (x, fx, n.pts = 256, ret = FALSE)
{

  if (class(fx) == "function")
    fx = fx(x)
  n.x = length(x)
  if (n.x != length(fx))
    stop("Unequal input vector lengths")
  if (n.pts < 64)
    n.pts = 64
  ap = stats::approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] +
                    ap$y[2 * (1:n.pts) + 1])/3
  invisible(list(value = sum(integral), cdf = list(x = ap$x[2 *
                                                              (1:n.pts)], y = cumsum(integral))))
}

# Author: Gery Geenens, extended by Kate Doubleday
#' Get KDE using Geenens et. al 2014 and 2018's methods: 1) transformation with probit or log function,
#' 2) local likelihood estimation
#' 3) nearest-neighbor bandwdith selection
#' @param x A vector of samples
#' @param xmin Minimum allowable x value (e.g., 0)
#' @param xmax Maximum allowable x value for probit transformation, or NA for log transformation (non-negative)
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @param max_scaler For log transforms, estimation is made over a range of width (max(x)-xmin)*max_scaler
#' @param weight One of 'LSCV', 'WLSCV1', 'WLSCV2'
#' @param n.res length of resulting estimation vector
#' @return A list of the estimate evaluation points, density, and cumulative distribution
probtranskde <- function(x, xmax, xmin=0, scale=0.9999, zero_offset=0.0001, max_scaler=2,
                         weight="WLSCV2",n.res=500) {

  if (!((is.na(xmax)|is.numeric(xmax)) & is.numeric(xmin))) stop('Bad input. Requires either xmin (xmax=NA) for log transform or xmin and xmax for probit transform')
  if (!(is.na(xmax)) & xmax <= xmin) stop("Bad input. Xmax must be greater than xmin.")

  xseq <- get_output_seq(x, xmin, xmax, n.res, max_scaler)

  n<-length(x)

  Ss <- transform_to_real_line(x, xmin, xmax, scale, zero_offset)

  # Get bandwidth from ks (kernel smoothing) package
  h0S<-ks::hpi(Ss)

  # Do KDE on the transformed space
  fhatSbasSs<-ks::kde(Ss,h=h0S,eval.points=Ss)$estimate

  a0<-cbind(seq(0.15,2,by=0.01),0)
  if(weight=="LSCV") {
    kLSCVseq<-locfit::lscvplot(~Ss,alpha=a0,deg=2,kern="gauss",renorm=TRUE,maxk=512)$value
    kLSCV<-mean(a0[which(kLSCVseq==min(kLSCVseq,na.rm=TRUE)),1])
    fit<-locfit::locfit(~lp(Ss,nn=kLSCV),deg=2,kern='gauss',renorm=TRUE)
  }
  if(weight=="WLSCV1") {
    kLSCVW1seq<-locfit::lscvplot(~Ss,alpha=a0,deg=2,kern="gauss",renorm=TRUE,maxk=512,weights=sqrt(dnorm(Ss)/fhatSbasSs))$value
    kLSCVW1<-mean(a0[which(kLSCVW1seq==min(kLSCVW1seq,na.rm=TRUE)),1])
    fit<-locfit::locfit(~lp(Ss,nn=kLSCVW1),deg=2,kern='gauss',renorm=TRUE)
  }
  if(weight=="WLSCV2") {
    kLSCVW2seq<-locfit::lscvplot(~Ss,alpha=a0,deg=2,kern="gauss",renorm=TRUE,maxk=512,weights=(dnorm(Ss)/fhatSbasSs))$value
    kLSCVW2<-mean(a0[which(kLSCVW2seq==min(kLSCVW2seq,na.rm=TRUE)),1])
    fit<-locfit::locfit(~lp(Ss,nn=kLSCVW2),deg=2,kern='gauss',renorm=TRUE)
  }

  # Predict is a generic function
  if (is.na(xmax)){
    ftildeX <- predict(fit, log(xseq))/xseq
  } else {
    ftildeX <- predict(fit,stats::qnorm(xseq))/stats::dnorm(stats::qnorm(xseq))
  }

  ftildeX<-ftildeX/myintegral(xseq,ftildeX)$value

  # Shift/scale back to the full spectrum
  res <- scale_full(xseq, ftildeX, xmin, xmax, scale)

  return(list(x=res$x, d=res$d, u=pracma::cumtrapz(res$x, res$d)))

}

get_output_seq <- function(x, xmin, xmax, n.res, scaler) {
  # On (0, inf) for log function
  if (is.na(xmax)) {
    return(seq(1/(n.res+1),n.res/(n.res+1),length=n.res)*(max(x)-xmin)*scaler)
  } else { # On (0,1) for probit function
    return(seq(1/(n.res+1),n.res/(n.res+1),length=n.res))
  }
}

#' Transform values from bounded to unbounded support through probit or log transform
#' Author: Gery Geenens and Kate Doubleday
#' @param x A vector of samples
#' @param xmin Minimum allowable x value (e.g., 0)
#' @param xmax Maximum allowable x value.
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @return The vector of transformed samples
transform_to_real_line <- function(x, xmin, xmax, scale, zero_offset) {
  # If need be, scale input vector to (0,1)
  Xs <- scale_01(x, xmin, xmax, scale, zero_offset)

  if (is.na(xmax)){
    Ss <- log(Xs)
  } else {  # Probit transformation, truncating weird extreme outliers
    Ss<-stats::qnorm(Xs)
  }

  if (any(is.infinite(Ss)|is.na(Ss))) stop('Transformation returned infinite or NA values. Check boundaries.')
  return(Ss)
}

#' Scale from [xmin, xmax] to (0,1). Values exactly at 0 are shifted up by an offset to stabilize the boundary.
#' All other values are scaled to avoid the boundary at 1.
#' @param x A vector of samples
#' @param xmin Minimum allowable x value (e.g., 0)
#' @param xmax Maximum allowable x value. If NA, function returns the input values.
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @return The vector of samples rescaled to [0,1]
scale_01 <- function(x, xmin, xmax, scale, zero_offset) {
  if (scale <=0 | scale >1) stop('Scale factor should be (0,1].')
  if (is.na(xmax)) {
    s <- x - xmin
    s[s==0] <- zero_offset
    return(s)
  } else {
    s <- (x - xmin)/xmax * scale
    s[s==0] <- zero_offset
    return(s)
  }
}

#' Scale from (0,1) to [xmin, xmax]. The scaling value is recovered in this back-transformation; the slight offset at 0 is not.
#' @param x01 A vector of points in (0,1)
#' @param d01 Associated density
#' @param xmin Minimum allowable x value (e.g., 0)
#' @param xmax Maximum allowable x value. If NA, function returns the input values.
#' @param scale Scaling factor (0,1]
#' @return The vector of samples rescaled to [0,1]
scale_full <- function(x01, d01, xmin, xmax, scale) {
  if (scale <=0 | scale >1) stop('Scale factor should be (0,1].')
  if (is.na(xmax)) {
    x <- x01 + xmin
    d <- d01
  } else {
    x <- x01/scale*xmax + xmin
    d <- d01/xmax*scale
  }
  return(list(x=x, d=d))
}

