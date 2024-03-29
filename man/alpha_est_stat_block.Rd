% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha_est_stat.R
\name{alpha_est_stat_block}
\alias{alpha_est_stat_block}
\title{Estimate the latent positions under the stationary assumption with block structure imposed}
\usage{
alpha_est_stat_block(A0, gamma_hat, d, K, lb, ub, sim = NULL)
}
\arguments{
\item{A0}{(n-by-n matrix) Initial adj. matrix}

\item{gamma_hat}{The estimated autoregressive parameter}

\item{d}{dimension of the latent position}

\item{K}{number of blocks}

\item{lb}{= -10 lower bound for \eqn{\hat{\Omega}_{ij}}}

\item{ub}{= 10 upper bound for \eqn{\hat{\Omega}_{ij}}}

\item{sim}{#' @return A list containing the following:
\describe{
\item{`alpha_est`}{(n-by-d matrix) Estimated latent positions}
\item{`Omega_hat`}{(n-by-n matrix) Estimated \eqn{\Omega} matrix}
}}
}
\description{
Estimate the latent positions under the stationary assumption with block structure imposed
}
