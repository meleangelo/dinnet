#' Estimate \eqn{\hat{\Omega}} (the matrix of fixed effects)
#'
#' @param d dimension of the latent position
#' @param A list of (n-by-n) sparse network matrices (one per time period)
#' @param gamma_hat estimated value of gamma from CMLE step
#'
#' @return A matrix
#' \describe{
#' \item{`Omega_hat`}{(n-by-n matrix) Estimated \eqn{\Omega} matrix}
#' }
#'
#' @export
#'
omega_est <- function(A, gamma_hat, d) {

  TT <- length(A)
  n = dim(A[[1]])[1]
  Omega_hat <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  for (t in 2:TT){
    tsvd <- irlba(A[[t]], d)
    Ipq <- dinnet::getIpq(A[[t]], d)
    Xhat <- tsvd$u %*% sqrt(diag(tsvd$d, ncol = d, nrow = d))
    rm(tsvd)
    gc()
    P_hat <- Xhat %*% Ipq %*% t(Xhat)
    P_hat <- grdpg::BFcheck(P_hat,l = 0.0001, u = 0.9999)
    rm(Xhat)
    gc()
    update <- ( log(P_hat / (1 - P_hat)) - gamma_hat * A[[t-1]] )/(TT-1)
    rm(P_hat)
    gc()
    Omega_hat <- Omega_hat + update
    rm(update)
    gc()

  }

  return(Omega_hat)
}
