% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{CVAR.fc_vine}
\alias{CVAR.fc_vine}
\title{Calculate VaR and CVaR from sampled data. CVaR calculation is done directly from the samples, rather than estimated from a fitted distribution.}
\usage{
\method{CVAR}{fc_vine}(x, samples = NA, epsilon = c(0.05, 0.95))
}
\arguments{
\item{x}{fc_vine object}

\item{samples}{(optional) previously obtained samples to use instead of new sampling, e.g. for coordination with quantiles calculation}

\item{epsilon}{Probability levels for lower/upper tail VaR/CVaR calculations, defaults to c(0.05, 0.95)}
}
\value{
list of var, cvar
}
\description{
Calculate VaR and CVaR from sampled data. CVaR calculation is done directly from the samples, rather than estimated from a fitted distribution.
}
