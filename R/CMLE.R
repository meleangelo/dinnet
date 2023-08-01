#' Estimate \eqn{\gamma} by CMLE
#'
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param gamma_init initial value of gamma
#' @param lb lower bound for gamma
#' @param ub upper bound for gamma
#'
#' @return \eqn{\hat{gamma}}
#'
#' @export
#'

CMLE_est <- function(A, gamma_init = NULL, xtol = 1.0e-5, lb = -1, ub = 1) {

  if (is.null(gamma_init)) gamma_init <- 0.5

  opts <- list("algorithm" =  "NLOPT_LD_LBFGS", #"NLOPT_LN_COBYLA",
               "xtol_rel" = xtol)

  opt_fit <- nloptr::nloptr(x0 = gamma_init,
                            eval_f = loglike_CMLE,
                            eval_grad_f = deriv_loglike_CMLE,
                            lb = lb,
                            ub = ub,
                            opts = opts,
                            A = A)

  return(opt_fit$solution)
}

#' Log likelihood function for the conditional maximum likelihood estimation
#'
#' @param gamma The autoregressive parameter
#' @param A (list of length T, containing n-by-n sparse adjacency matrices) time series of the adjacency matrices
#'
#' @return log-likelihood (times -1)
#'

loglike_CMLE <- function(gamma, A) {

  # read in dimensions
  n <- dim(A[[1]])[1]
  TT <- length(A)
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Initialize the log-likelihood
  loglik <- 0

  for (t in 2:TT2) {
    for (s in (t + 1):TT1) {

      # compute the common factor in both numerator and denominator
      if (s - t <= 2) {
        phi <- gamma * (A[[t - 1]] - A[[s + 1]])
      } else {
        phi <- gamma * (A[[t - 1]] - A[[s + 1]] +
                        A[[t + 1]] - A[[s - 1]])
      }



      # Compute the likelihood contribution for the (s,t) pair
      used_links = Matrix::Matrix((A[[t]] + A[[s]] == 1), sparse=TRUE)
      #loglik_t <- 0.5 * sum(
      # (A[[t]] + A[[s]] == 1) * (
      #   (A[[t]] * phi) - log(1 + exp(phi))
      # )
      #) / (n^2)

      loglik_t <- 0.5 * sum( used_links * (A[[t]] * phi) ) / (n^2)

      # compute log(1+exp(phi)) in a smart way to avoid memory problems
      temp_matrix = Matrix::Matrix(phi[used_links ], sparse = TRUE)
      logsum = 0.5 *sum(log(1+exp(temp_matrix)))  / (n^2)

      loglik_t = loglik_t - logsum

      #phi_values = unique(as.vector(phi))
      #compute_sparse = 0
      #for (p in 1:length(phi_values)){
      #  compute_sparse = compute_sparse + 0.5 * sum(phi == phi_values[p])*log( 1 + exp(phi_values[p]))/(n^2)
      #}

      #loglik_t = loglik_t - compute_sparse

      loglik <- loglik + loglik_t
    }
  }
  return(-loglik)
}

#' Derivative of the Log likelihood function for the conditional
#' maximum likelihood estimation
#'
#' @importFrom grdpg sigmoid
#'
#' @inheritParams loglike_CMLE
#'
#' @return derivative of the log-likelihood (times -1)

deriv_loglike_CMLE <- function(gamma, A) {

  # read in dimensions
  n <- dim(A[[1]])[1]
  TT <- length(A)
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Intialize the derivative of the log-likelihood
  deriv_loglik <- 0

  for (t in 2:TT2) {
    for (s in (t + 1):TT1) {

      if (s - t <= 2) {
        deriv_phi <- A[[t - 1]] - A[[s + 1]]

      } else {
        deriv_phi <- A[[t - 1]] - A[[s + 1]] + A[[t + 1]] - A[[s - 1]]
      }
      phi <- gamma * deriv_phi

      # Compute the derivative of the likelihood contribution for the (s,t) pair
      # deriv_loglik_t <- 0.5 * sum(
      #   (A[[t]] + A[[s]] == 1) * (A[[t]] - sigmoid(phi)) * deriv_phi
      # ) / (n^2)

      used_links = Matrix::Matrix((A[[t]] + A[[s]] == 1), sparse=TRUE)
      deriv_loglik_t = 0.5 * sum(used_links * A[[t]] * deriv_phi) / (n^2) # first part
      phi_sparse = Matrix::Matrix(phi[used_links ], sparse = TRUE)
      deriv_phi_sparse = Matrix::Matrix(deriv_phi[used_links ], sparse = TRUE)
      second_part = 0.5 *sum( sigmoid(phi_sparse) * deriv_phi_sparse )  / (n^2)
      deriv_loglik_t = deriv_loglik_t - second_part


      deriv_loglik <- deriv_loglik + deriv_loglik_t
    }
  }
  return(-deriv_loglik)
}


