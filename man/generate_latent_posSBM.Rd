% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{generate_latent_posSBM}
\alias{generate_latent_posSBM}
\title{Function to generate K latent positions from a stochastic block model}
\usage{
generate_latent_posSBM(latent, d, block_size)
}
\arguments{
\item{latent}{(d-by-K Matrix) latent position for each block}

\item{d}{(Int) Dimension of the latent position}

\item{block_size}{(Int vector of length K)}
}
\value{
(n-by-d matrix) latent position matrix (\eqn{\alpha} matrix)
  where n is number of nodes
}
\description{
The function generates latent position assuming \eqn{a_i} takes K values (in
the matrix latent), i.e. from a stochastic block model
}
\examples{
latent <- matrix(c(1, -2), 1, 2)
block_size <- c(10, 20)
latent_positions <- generate_latent_posSBM(latent, 1, block_size)

}
