#' Estimate \eqn{\gamma} by CMLE
#'
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param gamma_init initial value of gamma
#'
#' @return \eqn{\hat{gamma}}
#'
#' @export
#'

CMLE_est_old <- function(A, gamma_init = NULL, xtol = 1.0e-5) {

  if (is.null(gamma_init)) gamma_init <- 0.5

  opts <- list("algorithm" = "NLOPT_LD_LBFGS",
               "xtol_rel" = xtol)

  opt_fit <- nloptr::nloptr(x0 = gamma_init,
                            eval_f = loglike_CMLE_old,
                            eval_grad_f = deriv_loglike_CMLE_old,
                            lb = -1,
                            ub = 1,
                            opts = opts,
                            A = A)

  return(opt_fit$solution)
}

#' Log likelihood function for the conditional maximum likelihood estimation
#'
#' @param gamma The autoregressive parameter
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#'
#' @return log-likelihood (times -1)
#'

loglike_CMLE_old <- function(gamma, A) {

  # read in dimensions
  n <- dim(A)[1]
  TT <- dim(A)[3]
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Initialize the log-likelihood
  loglik <- 0

  for (t in 2:TT2) {
    for (s in (t + 1):TT1) {

      # compute the common factor in both numerator and denominator
      if (s - t <= 2) {
        phi <- gamma * (A[, , (t - 1)] - A[, , (s + 1)])
      } else {
        phi <- gamma * (A[, , (t - 1)] - A[, , (s + 1)] +
                          A[, , (t + 1)] - A[, , (s - 1)])
      }

      # Compute the likelihood contribution for the (s,t) pair
      loglik_t <- 0.5 * sum(
        (A[, , t] + A[, , s] == 1) * (
          (A[, , t] * phi) - log(1 + exp(phi))
        )
      ) / (n^2)

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

deriv_loglike_CMLE_old <- function(gamma, A) {

  # read in dimensions
  n <- dim(A)[1]
  TT <- dim(A)[3]
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Intialize the derivative of the log-likelihood
  deriv_loglik <- 0

  for (t in 2:TT2) {
    for (s in (t + 1):TT1) {

      if (s - t <= 2) {
        deriv_phi <- A[, , (t - 1)] - A[, , (s + 1)]

      } else {
        deriv_phi <- A[, , (t - 1)] - A[, , (s + 1)] + A[, , (t + 1)] - A[, , (s - 1)]
      }
      phi <- gamma * deriv_phi

      # Compute the derivative of the likelihood contribution for the (s,t) pair
      deriv_loglik_t <- 0.5 * sum(
        (A[, , t] + A[, , s] == 1) * (A[, , t] - sigmoid(phi)) * deriv_phi
      ) / (n^2)

      deriv_loglik <- deriv_loglik + deriv_loglik_t
    }
  }
  return(-deriv_loglik)
}


