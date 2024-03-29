% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMLE.R
\name{deriv_loglike_CMLE}
\alias{deriv_loglike_CMLE}
\title{Derivative of the Log likelihood function for the conditional
maximum likelihood estimation}
\usage{
deriv_loglike_CMLE(gamma, A)
}
\arguments{
\item{gamma}{The autoregressive parameter}

\item{A}{(n-by-n-by-T array) time series of the adjacency matrices}
}
\value{
derivative of the log-likelihood (times -1)
}
\description{
Derivative of the Log likelihood function for the conditional
maximum likelihood estimation
}
