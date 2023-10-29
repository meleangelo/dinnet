#' Estimate \eqn{\alpha} from \eqn{\hat{\Omega}}
#'
#' @param Omega_hat matrix of fixed effects
#' @param d dimension of the latent position
#'
#' @return A matrix
#' \describe{
#' \item{`alpha_hat`}{(n-by-d matrix) Estimated latent positions}
#' }
#'
#' @export
#'
alpha_est <- function(Omega_hat, d) {

  tsvd <- irlba(Omega_hat, d)

    if (length(tsvd$d) > 1){
      alpha_hat <- tsvd$u %*% sqrt(diag(tsvd$d, ncol = d, nrow = d))
    }
    else {
      alpha_hat <- tsvd$u * sqrt(tsvd$d)
    }
    return(alpha_hat)
}
