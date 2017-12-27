
#' @title Simulation for a Gamma Distribution with mean parametrization
#'
#' @description
#' Random generator for the Gamma distribution, using parameters \code{mean}
#' and \code{sigma}. This is an alternative to \code{rgamma} function.
#'
#' @details
#' The Gamma distribution with parameters \code{mean = mu} and \code{sigma = s} is
#' defined as
#' \deqn{f(x) = (s/mu)^s /Gamma(s) x^(s-1) e^-(s*x/mu)}
#' for x >= 0, mu > 0 and s > 0.
#'
#' The mean and variance of this parametrization are \eqn{E(x) = mu} and
#' \eqn{Var(X) = mu^2/sigma}
#' @param n Number of observations.
#' @param mu Mean parameter of the Gamma distribution.
#' @param sigma Scale parameter of the Gamma distribution.
#'
#' @return
#' A vector of length \code{n}, representing the realization of the Gamma
#' distribution with parameters \code{mean} and \code{sigma}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' rgamma_mu(100, 5, 5)
#'
#' @importFrom stats rgamma
#' @export
rgamma_mu <- function (n, mu, sigma) {
  rgamma(n, shape  = sigma, rate = sigma / mu)
}
