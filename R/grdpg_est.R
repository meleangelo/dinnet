#' Function for estimation of P with ASE period by period
#'
#'
#' @section Warning: random seeds changed
#'
#' @param A_t (n-by-n) Adjacency matrix
#' @param d dimension of the latent position
#' @param sim The random seed
#'
#' @return estimated P
#'
#' @export

estimate_ASE <- function(A_t, d, sim = NULL, bfcheck = TRUE, l = 0.0001, u = 0.9999){
  # Set the random seeds for the DGP
  if (!is.null(sim)) set.seed(sim)
  # estimation of P_t
  ASE <- grdpg::SpectralEmbedding(A_t, d)
  Ipq <- grdpg::getIpq(A_t, d)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d))
  Phat <- Xhat %*% Ipq %*% t(Xhat)

  if (bfcheck){
    return(grdpg::BFcheck(Phat, l = l, u = u))
  }  else {
    return(Phat)
  }
}

#' Function for estimation of X with ASE
#'
#'
#' @section Warning: random seeds changed
#'
#' @param A_t (n-by-n) Adjacency matrix
#' @param d dimension of the latent position
#' @param sim The random seed
#'
#' @return estimated X_hat
#'
#' @export
estimateX_ASE <- function(A, d, sim = NULL){
  # Set the random seeds for the DGP
  if (!is.null(sim)) set.seed(sim)
  # estimation of P_t
  ASE <- grdpg::SpectralEmbedding(A, d)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d) )
  return(Xhat)
}

#' Function for estimation of P with GRDPG
#'
#'
#' @section Warning: random seeds changed
#'
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param d dimension of the latent position
#' @param sim The random seed
#'
#' @return estimated P array for all time periods
#'
#' @export

estimate_GRDPG <- function(A, d, sim = NULL, bfcheck = TRUE, l = 0.0001, u = 0.9999){

  # read in dimensions
  n <- dim(A[[1]])[1]
  TT <- length(A)
  Phat <- vector("list", length = TT)
  for (t in 1:TT) {
    Phat[[t]] <- estimate_ASE(A[[t]], d, sim, bfcheck, l, u)
  }
  return(Phat)
}
