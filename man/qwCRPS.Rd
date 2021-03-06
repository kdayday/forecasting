% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_forecast.R
\name{qwCRPS}
\alias{qwCRPS}
\title{Get the quantile-weighted CRPS}
\usage{
qwCRPS(ts, tel, weighting = "none", quantiles = seq(0.01, 0.99, by =
  0.01), qs = NA)
}
\arguments{
\item{ts}{A ts_forecast object}

\item{tel}{A vector of the telemetry values}

\item{weighting}{One of "none" (default), "tails", "right", "left", "center"}

\item{quantiles}{(optional) Sequence of quantiles to integrate over}

\item{qs}{(optional) A list of quantile scores corresponding to the quantiles}
}
\description{
Calculate CRPS (continuous ranked probability score) for the forecast
estimated by trapezoidal integration of its quantile decomposition Quantile
weighting functions as suggested in Gneiting & Ranjan 2012
}
