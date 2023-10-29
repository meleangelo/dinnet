#' Get `Ipq` Matrix for GRDPG
#'
#' Construct `Ipq` matrix for Generalized Random Dot Product Graph.
#'
#' @import RSpectra
#'
#' @param A A square matrix.
#' @param d Embeded dimension.
#' @param method Method used to compute eigenvalues. Coule be `RSpectra` (default) or `eigen`.
#'
#' @return A `d` by `d` diagonal matrix where `d` is the embeded dimension with 1 and -1 on the diagonal.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @references Rubin-Delanchy, P., Priebe, C. E., Tang, M., & Cape, J. (2017). A statistical interpretation of spectral embedding: the generalised random dot product graph. \emph{arXiv preprint \href{https://arxiv.org/abs/1709.05506}{arXiv:1709.05506}}.
#'
#' @seealso \code{\link{grdpg}}, \code{\link{eigs}}, \code{\link{eigen}}
#'
#' @export


getIpq <- function(A, d) {

  if (d == 1) {
    Ipq <- matrix(1)
  } else {
    cols <- ncol(A)
    temp1 <- eigs_sym(A, d, 'LA')
    s1 <- temp1$values
    temp2 <- eigs_sym(A, d, 'SA')
    s2 <- temp2$values
    tempdat <- data.frame(raw = c(s1,s2)) %>%
      mutate(sign = ifelse(raw>0, 1, -1), s = abs(raw)) %>%
      arrange(desc(s))
    Ipq <- diag(tempdat$sign[1:d])
  }

  return(Ipq)
}
