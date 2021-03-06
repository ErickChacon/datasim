% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effects-multi.R
\name{mfe}
\alias{mfe}
\title{Multivariate Fixed Effects}
\usage{
mfe(x, beta)
}
\arguments{
\item{x}{A replicated continuous covariate of length \code{nq} for which the
fixed effects are evaluated. It is a \eqn{q}-times replicated vector of \eqn{x*}.}

\item{beta}{A vector of length \code{q}, where each element is a regression
coefficient associated to each response variable.}
}
\value{
A \eqn{nq}-length numeric vector of the evaluated multivariate effects.
}
\description{
\code{mfe} evaluates the multivariate fixed effects for a continuous replicated
covariate \code{x}.
}
\details{
Considering \eqn{x*} the continuous covariate of length \eqn{n}, the returning
multivariate effect is
\deqn{vec(x*b),}
where \eqn{b} is a \eqn{q}-length row vector of regression coefficients and
\eqn{vec(.)} represents a vectorization by columns of the provided matrix.
}
\examples{
# Different effects for each response.
mfe(x = rep(rnorm(5), 3), beta = c(0.1, 0, 1))

# Same effects for each response.
mfe(x = rep(rnorm(5), 3), beta = rep(10, 3))

}
\author{
Erick A. Chacón-Montalván
}
