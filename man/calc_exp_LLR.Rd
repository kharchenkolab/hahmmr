% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflows.R
\name{calc_exp_LLR}
\alias{calc_exp_LLR}
\title{Calculate LLR for an expression HMM}
\usage{
calc_exp_LLR(
  Y_obs,
  lambda_ref,
  d,
  phi_mle,
  mu = NULL,
  sig = NULL,
  alpha = NULL,
  beta = NULL
)
}
\arguments{
\item{Y_obs}{numeric vector Gene expression counts}

\item{lambda_ref}{numeric vector Reference expression levels}

\item{d}{numeric vector Total library size}

\item{phi_mle}{numeric MLE of expression fold change phi (alternative hypothesis)}

\item{mu}{numeric Mean parameter for the PLN expression model}

\item{sig}{numeric Dispersion parameter for the PLN expression model}

\item{alpha}{numeric Hyperparameter for the gamma poisson model (not used)}

\item{beta}{numeric Hyperparameter for the gamma poisson model (not used)}
}
\value{
numeric Log-likelihood ratio
}
\description{
Calculate LLR for an expression HMM
}
\keyword{internal}
