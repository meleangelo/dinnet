#' Function to generate K latent positions from a stochastic block model
#'
#' The function generates latent position assuming \eqn{a_i} takes K values (in
#' the matrix latent), i.e. from a stochastic block model
#'
#' @param latent (d-by-K Matrix) latent position for each block
#' @param d (Int) Dimension of the latent position
#' @param block_size (Int vector of length K)
#'
#' @return (n-by-d matrix) latent position matrix (\eqn{\alpha} matrix)
#'   where n is number of nodes
#'
#' @examples
#' latent <- matrix(c(1, -2), 1, 2)
#' block_size <- c(10, 20)
#' latent_positions <- generate_latent_posSBM(latent, 1, block_size)
#'
#' @export
#'
generate_latent_posSBM <- function(latent, d, block_size){

  # Check inputs
  if (nrow(latent) != d) {
    stop("The number of rows in `latent` should equal to `d`.")
  }
  if (ncol(latent) != length(block_size)) {
    stop("The number of blocks (K) in 'block_size' and 'latent' should match.")
  }

  # read number of blocks
  K <- length(block_size)

  # Compute the latent position (n-by-d matrix)
  latent_positions <- t(latent[, rep(seq.int(1, K), block_size), drop = FALSE])
  return(latent_positions)
}

#' Function to generate the network time series
#'
#' The function generates the (stationary) time series of the networks
#' without covariates \eqn{x_{ij}}
#'
#'
#' @importFrom grdpg sigmoid
#' @import Matrix
#'
#' @section Warning: random seeds changed
#'
#' @param latent_positions (n-by-d Matrix) The alpha matrix
#' @param gamma (\eqn{|\gamma| < 1}) The autoregressive parameter
#' @param time_periods number of time periods T
#' @param sim The random seed
#'
#'
#' @return A list containing the following:
#' \describe{
#' \item{`A0`}{(n-by-n matrix) Initial adjacency matrix \eqn{A_0}}
#' \item{`P0`}{(n-by-n matrix) Initial transition probability matrix \eqn{P_0}}
#' \item{`A`}{(n-by-n-by-T array) time series of the adjacency matrices}
#' \item{`P`}{(n-by-n-by-T array) time series of the transition probability matrix}
#' \item{`alpha`}{(n-by-n matrix) the product of latent positions \eqn{a_i^T a_j}}
#' }
#'
#' @export

generate_data <- function(latent_positions, gamma, time_periods, sim = NULL){

  # Set the random seeds for the DGP
  if (!is.null(sim)) set.seed(sim)

  # read number of nodes
  n <- nrow(latent_positions)

  # generate fixed effects matrix
  alpha <- latent_positions %*% t(latent_positions)

  # generate initial matrix P_0 (stationary distribution)
  P0 <- sigmoid(alpha) / (1 - sigmoid(gamma + alpha) + sigmoid(alpha))

  # generate initial adjacency A_0
  A0 <- sym_mat_0(
    1 * (matrix(runif(n^2), ncol = n, nrow = n) < P0)
  )

  # initialize P and A
  P <- vector("list", length = time_periods)
  A <- vector("list", length = time_periods)

  P[[1]] <- sigmoid(gamma * A0 + alpha)
  A[[1]] <- Matrix::Matrix( sym_mat_0(1* (matrix(runif(n ^ 2), ncol = n, nrow = n) < P[[1]] )), sparse = TRUE)



  for (t in 2:time_periods) {
    t1 <- t - 1
    P[[t]] <- sigmoid(gamma * A[[t1]] + alpha)
    A[[t]] <- Matrix::Matrix( sym_mat_0(1* (matrix(runif(n ^ 2), ncol = n, nrow = n) < P[[t]] )), sparse = TRUE)

  }

  data <- list(A0, P0, A, P, alpha)
  names(data) <- c("A0", "P0", "A", "P", "alpha")
  return(data)
}


