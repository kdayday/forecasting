#' A structure for transforming back and forth to uniform domain
#'
#' @param x A vector of samples
#' @param method Method for estimating the marginal distribution. One of 'empirical', 'kde1d', 'geenens', 'kernsmooth'
#' @param ... Optional arguments to the distribution estimator. 'xmax' can be used for 'kde1d' and 'geenens'
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
  pdf <- kde1d::kde1d(x, ...)
  cdf <- pracma::cumtrapz(pdf$grid_points, pdf$values)
  return(list(x=pdf$grid_points, d=as.numeric(pdf$values), u=as.numeric(cdf)))
}

#' Get probability using empirical CDF, with a default minimum x value of 0
#'
#' @param x A vector of samples
#' @param xmax maximum
#' @return A list of the evaluation points, density, and cumulative distribution
probempirical <- function(x, xmax=max(ceiling(x))) {

  # Pick the largest maximum value, in case actual data exceeds rating
  xmax = max(xmax, max(ceiling(x)))
  if (xmax > max(x)) {
    xseq <- c(0, sort(x), xmax)
  } else {
      xseq <- c(0, sort(x))
  }

  # Aggregate to combine duplicate values
  agg <- stats::aggregate(xseq, list(value = xseq), length)
  xout <- as.numeric(unlist(agg['value']))
  u <- as.numeric(unlist((cumsum(agg['x'])-1)/(length(xseq)-1)))

  return(list(x=xout, u=u, d=NaN))
}

# ----------------------------------------------------------------------------
# Functions for bounded support transformations with nearest neighbors bandwidths

# Author: Gery Geenens
myintegral<-function (x, fx, n.pts = 256)
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
  invisible(list(value = sum(integral), cdf = list(x = ap$x[2 * (1:n.pts)], y = cumsum(integral))))
}

# Author: Gery Geenens, extended by Kate Doubleday
#' Get KDE using Geenens et. al 2014 and 2018's methods: 1) transformation with probit or log function,
#' 2) local likelihood estimation
#' 3) nearest-neighbor bandwdith selection
#' @param x A vector of samples
#' @param xmax Maximum allowable x value for probit transformation, or NaN for log transformation (non-negative)
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @param max_scaler For log transforms, estimation is made over the range max(x)*max_scaler
#' @param weight One of 'LSCV', 'WLSCV1', 'WLSCV2'
#' @param n.res length of resulting estimation vector
#' @return A list of the estimate evaluation points, density, and cumulative distribution
probtranskde <- function(x, xmax, scale=0.9999, zero_offset=0.0001, max_scaler=2,
                         weight="WLSCV2",n.res=500) {

  if (!((is.nan(xmax)|is.numeric(xmax)))) stop('Bad input. Requires either xmax=NaN for log transform or numeric xmax for probit transform')
  if (!(is.nan(xmax)) & xmax <= 0) stop("Bad input. Xmax must be non-negative.")

  xseq <- get_output_seq(x, xmax, n.res, max_scaler)

  n<-length(x)

  Ss <- transform_to_real_line(x, xmax, scale, zero_offset)

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
  if (is.nan(xmax)){
    ftildeX <- predict(fit, log(xseq))/xseq
  } else {
    ftildeX <- predict(fit,stats::qnorm(xseq))/stats::dnorm(stats::qnorm(xseq))
  }

  ftildeX<-ftildeX/myintegral(xseq,ftildeX)$value

  # Shift/scale back to the full spectrum
  return(scale_full(xseq, ftildeX, xmax))
}

get_output_seq <- function(x, xmax, n.res, scaler) {
  # On (0, inf) for log function
  if (is.nan(xmax)) {
    return(seq(1/(n.res+1),n.res/(n.res+1),length=n.res)*max(x)*scaler)
  } else { # On (0,1) for probit function
    return(seq(1/(n.res+1),n.res/(n.res+1),length=n.res))
  }
}

#' Transform values from bounded to unbounded support through probit or log transform. Assumes lower boundary is at 0.
#' Author: Gery Geenens and Kate Doubleday
#' @param x A vector of samples
#' @param xmax Maximum allowable x value.
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @return The vector of transformed samples
transform_to_real_line <- function(x, xmax, scale, zero_offset) {
  # If need be, scale input vector to (0,1)
  Xs <- scale_01(x, xmax, scale, zero_offset)

  if (is.nan(xmax)){
    Ss <- log(Xs)
  } else {  # Probit transformation, truncating weird extreme outliers
    Ss<-stats::qnorm(Xs)
  }

  if (any(is.infinite(Ss)|is.na(Ss))) stop('Transformation returned infinite or NA/NaN values. Check boundaries.')
  return(Ss)
}

#' Scale from [0, xmax] to (0,1). Values exactly at 0 are shifted up by an offset to stabilize the boundary.
#' All other values are scaled to avoid the boundary at 1.
#' @param x A vector of samples
#' @param xmax Maximum allowable x value. If NA, function returns the input values.
#' @param scale Scaling factor (0,1] to move maximum values off boundary at 1
#' @param zero_offset Amount to shift minimum values off the boundary at 0 in the (0,1) domain
#' @return The vector of samples rescaled to [0,1]
scale_01 <- function(x, xmax, scale, zero_offset) {
  if (scale <=0 | scale >1) stop('Scale factor should be (0,1].')
  if (is.nan(xmax)) {s <- x} else {s <- x/xmax * scale}
  s[s==0] <- zero_offset
  return(s)
}

#' Scale from (0,1) to [0, xmax]. The slight offset at 0 is not recovered here.
#' @param x01 A vector of points in (0,1)
#' @param d01 Associated density
#' @param xmax Maximum allowable x value. If NA, function returns the input values.
#' @return The vector of samples rescaled to [0,1]
scale_full <- function(x01, d01, xmax) {
  if (x01[1] > 0) { # Add the hard boundary at 0 back in if needed
    x01 <- c(0, x01)
    d01 <- c(0, d01)
  }
  if (is.nan(xmax)) {
      x <- x01
      d <- d01
  } else {
    scaler <- xmax/max(x01)
    x <- x01*scaler
    d <- d01/scaler
  }
  u <- pracma::cumtrapz(x, d)
  return(list(x=x, d=d/max(u), u=u/max(u)))
}

