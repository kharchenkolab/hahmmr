% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflows.R
\name{fit_ref_sse}
\alias{fit_ref_sse}
\title{Fit a reference profile from multiple references using constrained least square}
\usage{
fit_ref_sse(Y_obs, lambdas_ref, gtf, min_lambda = 2e-06, verbose = FALSE)
}
\arguments{
\item{Y_obs}{vector}

\item{lambdas_ref}{matrix}

\item{gtf}{dataframe}
}
\value{
fitted expression profile
}
\description{
Fit a reference profile from multiple references using constrained least square
}
\keyword{internal}
