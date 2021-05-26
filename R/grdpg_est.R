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

estimate_ASE <- function(A_t, d, sim){
  set.seed(sim)
  # estimation of P_t
  ASE <- grdpg::SpectralEmbedding(A_t, d, work  = 100)
  Ipq <- grdpg::getIpq(A_t, d)
  Xhat <- ASE$X %*% sqrt(diag(ASE$D, ncol = d, nrow = d))
  Phat <- Xhat %*% Ipq %*% t(Xhat)
  return(Phat)
}


#' Function for estimation of P with GRDPG
#'
#'
#' @section Warning: random seeds changed
#'
#' @inheritParams estimate_ASE
#'
#' @return estimated P array for all time periods
#'
#' @export

estimate_GRDPG <- function(A, d, sim){

  set.seed(sim)

  n <- dim(A)[1]
  TT <- dim(A)[3]
  Phat <- array(NA, dim = c(n,n,TT))
  for (t in 1:TT){
    Phat[, , t] <- estimate_ASE(A[, , t], d, sim)
  }
  return(Phat)
}
