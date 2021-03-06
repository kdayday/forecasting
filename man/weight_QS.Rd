% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_forecast.R
\name{weight_QS}
\alias{weight_QS}
\title{Calculate a vector of weighted quantile scores}
\usage{
weight_QS(qs, quantiles, weighting = "none")
}
\arguments{
\item{qs}{A vector of quantile scores}

\item{quantiles}{A vector of the quantiles (percentiles) in [0,1]}

\item{weighting}{One of "none" (default), "tails", "right", "left", "center"}
}
\value{
A vector of the weighted scores
}
\description{
Can emphasize one or both tails or center
}
