% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bma.R
\name{get_lm}
\alias{get_lm}
\title{Get the linear model for a single ensemble member's 'b' coefficients
Telemetry within given tolerance of 1 are assumed to be clipped glm and lm
omit NA values by default}
\usage{
get_lm(fc, tel, form, B_transform, percent_clipping_threshold)
}
\arguments{
\item{fc}{Vector of training forecast data on [0,1]}

\item{tel}{Vector of telemetry data on [0,1]}

\item{form}{Formula in terms of x,y for logistic regression model}

\item{B_transform}{A function for pre-transforming the forecast data
(optional), e.g. function(X) return(log(x)+1)}

\item{percent_clipping_threshold}{Percentage threshold for determining if
clipping is occurring}
}
\value{
Summary of the lm model
}
\description{
Get the linear model for a single ensemble member's 'b' coefficients
Telemetry within given tolerance of 1 are assumed to be clipped glm and lm
omit NA values by default
}
