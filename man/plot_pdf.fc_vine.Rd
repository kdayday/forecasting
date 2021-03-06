% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{plot_pdf.fc_vine}
\alias{plot_pdf.fc_vine}
\title{Plot probabilistic forecast's estimated pdf (with kde)
Note that CVaR and VaR, while represented on the graph, are calculated directly from sampled data rather than estimated
from the kde results.}
\usage{
\method{plot_pdf}{fc_vine}(x, cvar = FALSE, epsilon = c(0.05, 0.95))
}
\arguments{
\item{x}{fc_vine object}
}
\description{
Plot probabilistic forecast's estimated pdf (with kde)
Note that CVaR and VaR, while represented on the graph, are calculated directly from sampled data rather than estimated
from the kde results.
}
