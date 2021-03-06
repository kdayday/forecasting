% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{fc_bma}
\alias{fc_bma}
\title{A BMA forecast}
\usage{
fc_bma(data.input, location, time, model, max_power, bma_distribution, ...)
}
\arguments{
\item{data.input}{A vector of ensemble members}

\item{location}{A string}

\item{time}{A lubridate time stamp}

\item{model}{A pre-fit BMA model from bma_ens_models}

\item{max_power}{Maximum power for normalizing forecast to [0,1]}

\item{bma_distribution}{One of "beta", "truncnorm" --> determines the type of
distribution used for each member's kernel dressing}

\item{...}{Additional parameters}
}
\value{
A probabilistic forecast object
}
\description{
Initialize a univariate probabilistic power forecast for a specific time
point using Bayesian model averaging (BMA)
}
\seealso{
Other forecast classes: \code{\link{fc_binned}},
  \code{\link{fc_emos}}, \code{\link{fc_empirical}}
}
\concept{forecast classes}
