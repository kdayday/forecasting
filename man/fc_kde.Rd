% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{fc_kde}
\alias{fc_kde}
\title{Initialize a univariate probabilistic power forecast for a specific time point using kernel density estimation. See kde_methods.R for more details}
\usage{
fc_kde(data.input, location, time, cdf.method = "geenens", ...)
}
\arguments{
\item{data.input}{A vector of ensemble members}

\item{location}{A string}

\item{time}{A lubridate time stamp}

\item{cdf.method}{KDE method selection, see kde_methods.R for details}

\item{...}{Additional parameters passed on the KDE method}
}
\value{
A probabilistic forecast object
}
\description{
Initialize a univariate probabilistic power forecast for a specific time point using kernel density estimation. See kde_methods.R for more details
}
