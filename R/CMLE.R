#' Log likelihood function for the conditional maximum likelihood estimation
#'
#' @param A (n-by-n-by-T array) time series of the adjacency matrices
#' @param gamma The autoregressive parameter
#'
#' @return log-likelihood
#'
#' @export

loglike_CMLE <- function(A, gamma) {

  # read in dimensions
  n <- dim(A)[1]
  TT <- dim(A)[3]
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Intialize the log-likelihood
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
  return(loglik)
}




#' Derivative of the Log likelihood function for the conditional
#' maximum likelihood estimation
#'
#' @importFrom grdpg sigmoid
#'
#' @inheritParams loglike_CMLE
#'
#' @return derivative of the log-likelihood

deriv_loglike_CMLE <- function(A, gamma) {

  # read in dimensions
  n <- dim(A)[1]
  TT <- dim(A)[3]
  TT2 <- TT - 2
  TT1 <- TT - 1

  # Intialize the log-likelihood
  deriv_loglik <- 0

  for (t in 2:TT2) {
    for (s in (t + 1):TT1) {

      if (s - t <= 2) {
        deriv_phi <- A[, , (t - 1)] - A[, , (s + 1)]

      } else {
        deriv_phi <- A[, , (t - 1)] - A[, , (s + 1)] + A[, , (t + 1)] - A[, , (s - 1)]
      }
      phi <- gamma * deriv_phi

      # Compute the likelihood contribution for the (s,t) pair
      deriv_loglik_t <- 0.5 * sum(
        (A[, , t] + A[, , s] == 1) * (
          (A[, , t] - sigmoid(phi)) * deriv_phi
        )
      ) / (n^2)

      deriv_loglik <- deriv_loglik + deriv_loglik_t
    }
  }
  return(deriv_loglik)
}
