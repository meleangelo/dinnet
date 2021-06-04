#' Estimate the latent positions under the stationary assumption using the initial adj. matrix
#'
#' @param A0 (n-by-n matrix) Initial adj. matrix
#' @param gamma_hat The estimated autoregressive parameter
#' @param d dimension of the latent position
#' @param lb = -10 lower bound for \eqn{\hat{\Omega}_{ij}}
#' @param ub = 10 upper bound for \eqn{\hat{\Omega}_{ij}}
#' @param sim
#'
#' @return A list containing the following:
#' \describe{
#' \item{`alpha_est`}{(n-by-d matrix) Estimated latent positions}
#' \item{`Omega_hat`}{(n-by-n matrix) Estimated \eqn{\Omega} matrix}
#' }
#'
#' @export
#'

alpha_est_stat <- function(A0, gamma_hat, d, lb, ub, sim = NULL) {

  P_hat <- estimate_ASE(A0, d, sim)

  Omega_hat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      sol <- uniroot(
        f = function(x){
          sigmoid(x) / (1 - sigmoid(gamma_hat + x) + sigmoid(x)) -  P_hat[i, j]
        }, interval = c(lb, ub)
      )$root
      Omega_hat[i, j] <- sol
    }
  }
  Omega_hat <- sym_mat(Omega_hat)

  svd_res <- grdpg::SpectralEmbedding(Omega_hat, d)
  if (length(svd_res$D) > 1) alpha_hat <- svd_res$X %*% sqrt(diag(svd_res$D))
  else alpha_hat <- svd_res$X * sqrt(svd_res$D)
  return(list(alpha_hat = alpha_hat,
              Omega_hat = Omega_hat))
}
