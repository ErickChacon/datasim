% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-model.R
\name{model_response}
\alias{model_response}
\title{Simulate response variable}
\usage{
model_response(model_frame, formula = attr(model_frame, "formula"),
  link_inv = replicate(length(formula), identity), generator = rnorm,
  responses = c("response"), seed = NULL)
}
\arguments{
\item{model_frame}{A \code{tibble} containing all the covariates to be used.
Usually this is obtained as the output of the \code{model_frame} function.}

\item{formula}{An optional list of formulas to simulate the parameters of the
response variables.}

\item{link_inv}{A list of function representing the inverse link function for the
parameters.}

\item{generator}{Function to generate the response variables given the parameters}

\item{responses}{character vector indicating the names of the response variables}

\item{seed}{Seed to be defined with function \code{set.seed} to obtain reproducible
results}
}
\value{
return.
}
\description{
\code{function} description.
}
\details{
details.
}
\examples{


f <- list(
  mean ~ I(5 + 0.5 * x1 + 0.1 * x2 + 0.7 * id1),
  sd ~ I(x1)
)
(model_fr <- model_frame(f, n = 10))
(data <- model_response(model_fr, link = list(identity, exp)))

beta0 <- c(-1, 1)
# Structure of the model
formula <- list(
  mean ~ fa(sex, beta = get("beta0")),
  sd ~ I(0)
)
(model_fr <- model_frame(formula, n = 10))
(data <- model_response(model_fr, link = list(identity, exp)))

formula <- list(
  mean ~ mfe(x, beta = 1:2),
  sd ~ mfe(x1, beta = 1:2)
)
(model_fr <- model_frame(formula, n = 10))
#(data <- model_response(model_fr, link_inv = list(identity, exp), responses = 1:2))
#(data <- sim_model(formula, link_inv = list(identity, exp), n = 10, responses = 1:2))

}
\author{
Erick A. Chacon-Montalvan
}
