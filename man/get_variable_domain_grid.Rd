% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{get_variable_domain_grid}
\alias{get_variable_domain_grid}
\title{Get a grid of evaluation points in the variable domain, based on the range of the marginal distribution estimation}
\usage{
get_variable_domain_grid(x, k)
}
\arguments{
\item{x}{A fc_vine object}

\item{k}{Integer or vector of number of samples to take along each dimension. If given an integer, all dimensions have the same number of samples.}
}
\value{
a grid of points size k^n
}
\description{
Get a grid of evaluation points in the variable domain, based on the range of the marginal distribution estimation
}
