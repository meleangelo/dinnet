#' Convert a matrix to a symmetric matrix
#'
#' Assign the upper triangular matrix to the lower triangle
#'
#'
#' @param A a square matrix
#'
#' @return `A_sym`
#'
sym_mat <- function(A) {
  A[lower.tri(A)] <- t(A)[lower.tri(t(A))]

  return(A)
}

#' Convert a matrix to a symmetric matrix
#'
#' Assign the upper triangular matrix to the lower triangle and let diagonal = 0
#'
#' @param A a square matrix
#'
#' @return `A_sym`
#'
sym_mat_0 <- function(A) {
  A[lower.tri(A)] <- t(A)[lower.tri(t(A))]
  diag(A) <- 0

  return(A)
}
