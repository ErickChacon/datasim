
#' @title Multivariate Fixed Effect
#'
#' @description
#' \code{mfe} compute the multivariate fixed effect. Generally used with
#' \code{msim_model}.
#'
#' @details
#' details.
#'
#' @param x A vector of length n for which the fixed effect will be evaluated.
#' @param beta A vector of length q, this is the fixed effect for each response
#' variable.
#'
#' @return A matrix of dimension n x q.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' mfe(x = rnorm(10), beta = c(0.1, 0, 1))
#'
#' @export
mfe <- function (x, beta) {
  as.numeric(matrix(x) %*% matrix(beta, nrow = 1))
}

