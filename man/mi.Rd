% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effects-multi.R
\name{mi}
\alias{mi}
\title{Multivariate Intercepts}
\usage{
mi(beta, size)
}
\arguments{
\item{beta}{A \eqn{q}-length numeric vector of intercepts
for each response variable.}

\item{size}{A numeric value \eqn{n} representing the number of units.}
}
\description{
\code{mi} handles the evaluation of multivariate intercepts.
}
\examples{
# Different effects for each response.
mi(c(1:3), size = 10)

}
\author{
Erick A. Chacón-Montalván
}
