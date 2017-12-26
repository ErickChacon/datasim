
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

#' @title Simulate a Multivariate Gaussian process
#'
#' @description
#' \code{mgp} Simulate a Multivariate spatial Gaussian process known as linear model
#' of coregionalization Y(h) = AS(h), where S(h) is a vector of q independent Gaussian
#' processes.
#'
#' @details
#' details.
#'
#' @param s1 First coordinate
#'
#' @param s2 Second coordinate
#'
#' @param cov.model A character or function indicating the covariance function that
#' Should be used to compute the variance-covariance matrix
#'
#' @param variance A qxq matrix of non-spatial covariance.
#' @param nugget A qxq diagonal matrix of non-spatial noise.
#' @param phi A q-length vector of decay parameters.
#' @param kappa A q-length vector of kappa parameters if Matern spatial correlation
#' function is used.
#'
#' @return A vector of the realization of the Gaussian Process
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' # Generate coordinates
#' N <- 100
#' s1 <- 2 * runif(N)
#' s2 <- 2 * runif(N)
#'
#' # Covariance parameters
#' q <- 2
#' var <- sqrt(diag(c(4, 4)))
#' A <- matrix(c(1, - 0.8, 0, 0.6), nrow = 2)
#' variance <- var %*% tcrossprod(A) %*% var
#' nugget <- diag(0, q)
#' phi <- rep(1 / 0.08, q)
#'
#' # Generate the multivariate Gaussian process
#' y <- mgp(s1, s2, "exponential", variance, nugget, phi)
#' y1 <- y[1:N]
#' y2 <- y[(N + 1):(2 * N)]
#'
#' # Check correlation
#' cor(y1, y2)
#' plot(y1, y2)
#'
#' # Visualize the spatial
#' plot(s1, s2, cex = y1, col = 2)
#' points(s1, s2, cex = y2, col = 3)
#'
#' @importFrom spBayes mkSpCov
#'
#' @export

mgp <- function (s1, s2, cov.model = NULL, variance = NULL, nugget = NULL, phi = NULL, kappa = NULL) {

  coords <- cbind(s1, s2)
  n <- nrow(coords)
  q <- nrow(variance)

  if (is.null(kappa)) {
    theta <- phi
  } else {
    theta <- c(phi, kappa)
  }

  varcov <- spBayes::mkSpCov(coords, variance, nugget, theta, cov.model)
  right <- chol(varcov)
  output <- as.numeric(crossprod(right, rnorm(n * q)))
  as.numeric(matrix(output, nrow = n, byrow = TRUE))

}
