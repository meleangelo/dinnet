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
  TT1 <- TT-1
  n <- dim(A[[1]])[1]

  # initialize matrix
  Omega_hat <- matrix(0, nrow = n, ncol = n)
  # loop to estimate matrix
  for (t in 2:TT){
    P_hat <- estimate_ASE(A[[t]], d, sim = 100)
    P_hat <- grdpg::BFcheck(P_hat)
    Omega_hat <- Omega_hat + ( log(P_hat / (1 - P_hat)) - gamma_hat * A[[t-1]] )/TT1
  }

  svd_res <- grdpg::SpectralEmbedding(Omega_hat, d)

  if (length(svd_res$D) > 1) alpha_hat <- svd_res$X %*% sqrt(diag(svd_res$D))
  else alpha_hat <- svd_res$X * sqrt(svd_res$D)
  return(list(alpha_hat = alpha_hat,
              Omega_hat = Omega_hat))
}

#' Brute-Force Check of Probability Matrix
#'
#' Check whether all entries of the estimated probabiltiy matrix are between 0 and 1.
#' Set the ones that not greater than 0 to be a lower bound and the ones that not less than 1 to be a upper bound.
#'
#' @param P Estimated probability matrix.
#' @param l Lower bound of the probability matrix, i.e. set all entries that not greater than 0 to `l`. 0.0001 by default.
#' @param u Upper bound of the probability matrix, i.e. set all entries that not less than 1 to be `u`. 0.9999 by default.
#'
#' @return A matrix of which all entries are between 0 and 1.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' P <- matrix(runif(25,-1,2), nrow = 5)
#' BFcheck(P)
#'
#' @export


BFcheck <- function(P, l = 0.0001, u = 0.9999) {
  if (!is.numeric(P)) {
    stop("Input need to be numeric.")
  }
  P[P<=0] <- l
  P[P>=1] <- u
  return(P)
}
