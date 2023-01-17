#' Estimate \eqn{\alpha} from \eqn{\hat{\Omega}}
#'
#' @inheritParams min_dist_omega
#' @param d dimension of the latent position
#' @param est_method estimation method of \eqn{\Omega} matrix. "avg" (direct estimation) or "opt" (via minimum distance)
#'
#' @return A list containing the following:
#' \describe{
#' \item{`alpha_est`}{(n-by-d matrix) Estimated latent positions}
#' \item{`Omega_hat`}{(n-by-n matrix) Estimated \eqn{\Omega} matrix}
#' }
#'
#' @export
#'
alpha_est <- function(A, gamma_hat, d, Omega_init = NULL, xtol = 1.0e-5) {
  # dimensions
  TT <- length(A)
  n <- dim(A[[1]])[1]

  # initialize matrix
  Omega_hat <- matrix(NA, nrow = n, ncol = n)
  # loop to estimate matrix
  for (t in 2:TT){
    P_hat <- estimate_ASE(A[[t]], d, sim = 100)
    Omega_hat <- Omega_hat + ( log(P_hat / (1 - P_hat)) - gamma_hat * A[[t-1]] )/TT
  }

  svd_res <- grdpg::SpectralEmbedding(Omega_hat, d)

  if (length(svd_res$D) > 1) alpha_hat <- svd_res$X %*% sqrt(diag(svd_res$D))
  else alpha_hat <- svd_res$X * sqrt(svd_res$D)
  return(list(alpha_hat = alpha_hat,
              Omega_hat = Omega_hat))
}

