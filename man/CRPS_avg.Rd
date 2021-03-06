% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_forecast.R
\name{CRPS_avg}
\alias{CRPS_avg}
\title{Get the average CRPS}
\usage{
CRPS_avg(ts, tel, normalize.by = 1, ...)
}
\arguments{
\item{ts}{A ts_forecast object}

\item{tel}{A vector of the telemetry values}

\item{normalize.by}{(optional) A normalization factor, either a scalar or
vector}

\item{...}{optional arguments to equalize_telemetry_forecast_length}
}
\description{
Get average estimated (continuous ranked probability score) for the
forecast; the score at each time point is estimated through numerical
integration. CRPS characterizes calibration and sharpness together
}
