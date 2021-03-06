% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_forecast.R
\name{get_transform_with_unique_xmin_max}
\alias{get_transform_with_unique_xmin_max}
\title{Subfunction for calc_transforms to cycle thorugh xmin and xmax if they are given uniquely for each dimension}
\usage{
get_transform_with_unique_xmin_max(idx, dat, cdf.method, ...)
}
\arguments{
\item{idx}{Column index of dat}

\item{dat}{training data matrix over all the dimensions}

\item{cdf.method}{marg_transform cdf.method}

\item{...}{Optional arguments to marg_transform, including potentially xmin or xmax in either scalar or vector form}
}
\description{
Subfunction for calc_transforms to cycle thorugh xmin and xmax if they are given uniquely for each dimension
}
