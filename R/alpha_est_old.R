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
alpha_est <- function(P_hat, A, gamma_hat, d, est_method = "avg", Omega_init = NULL, xtol = 1.0e-5) {

  if (est_method == "avg") {
    Omega_hat <- avg_omega(P_hat, A, gamma_hat)
  } else if (est_method == "opt") {
    Omega_hat <- min_dist_omega(P_hat, A, gamma_hat, Omega_init, xtol)
  } else {
    stop("Invalid method input.")
  }

  svd_res <- grdpg::SpectralEmbedding(Omega_hat, d)

  if (length(svd_res$D) > 1) alpha_hat <- svd_res$X %*% sqrt(diag(svd_res$D))
  else alpha_hat <- svd_res$X * sqrt(svd_res$D)
  return(list(alpha_hat = alpha_hat,
              Omega_hat = Omega_hat))
}

#' Estimation of \eqn{\hat{\Omega}} by averaging
#'
#' @param P_hat (n-by-n-by-T array) Estimated time series of the transition probability matrix
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param gamma_hat The estimated autoregressive parameter
#'
#' @return \eqn{\Omega} matrix
#'
#' @export
#'
#'
avg_omega <- function(P_hat, A, gamma_hat) {
  TT <- dim(A)[3]

  P_temp <- P_hat[, , 2:TT]
  A_temp <- A[, , 1:(TT - 1)]

  Omega <- apply(log(P_temp / (1 - P_temp)) - gamma_hat * A_temp, c(1,2), mean)

  return(Omega)
}


#' Minimize distance function to estimation \eqn{\hat{\Omega}}
#'
#' @param P_hat (n-by-n-by-T array) Estimated time series of the transition probability matrix
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param gamma_hat The estimated autoregressive parameter
#' @param Omega_init = NUll by default: initial value
#' @param xtol = 1.0e-5 by default: numerical accuracy
#'
#' @return \eqn{\Omega} matrix
#'
#' @export
min_dist_omega <- function(P_hat, A, gamma_hat, Omega_init = NULL, xtol = 1.0e-5) {
  # do the minimization entry-by-entry

  n <- dim(A)[1]
  TT <- dim(A)[3]

  if (is.null(Omega_init)) {
    Omega_init <- matrix(0, n, n) # Cannot verify convexity of the obj.
                                  # Think about how to generate the initial value
  }

  Omega <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {

      opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                   "xtol_rel" = xtol)
      opt_fit <- nloptr::nloptr(x0 = Omega_init[i, j],
                                eval_f = dist_ij,
                                eval_grad_f = deriv_dist_ij,
                                opts = opts,
                                p = P_hat[i, j, ],
                                a = A[i, j, ],
                                gamma = gamma_hat)
      Omega[i,j] = opt_fit$solution
    }
  }
  return(sym_mat(Omega))

}

#' entry-wise distance function
#'
#' @importFrom grdpg sigmoid
#'
#' @param omega scalar
#' @param p vector of length TT
#' @param a vector of length TT
#' @param gamma The autoregressive parameter
#'
#' @return distance_ij
#'
#'
dist_ij <- function(omega, p, a, gamma) {
  TT <- length(p)
  sum((p[2:TT] - sigmoid(gamma * a[1:(TT - 1)] + omega))^2) / (TT - 1)
}


#' Derivative of the entry-wise distance function
#'
#' @importFrom grdpg sigmoid
#'
#' @inheritParams dist_ij
#'
#' @return distance_ij
#'
#'
deriv_dist_ij <- function(omega, p, a, gamma) {
  TT <- length(p)
  x <- gamma * a[1:(TT - 1)] + omega
  p <- p[2:TT]
  2 * sum(
    (p - sigmoid(x)) * sigmoid(x) * (sigmoid(x) - 1)
  ) / (TT - 1)
}
