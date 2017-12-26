
#' @title Simulate databases based in models
#'
#' @description
#' \code{sim_model} simulate a database based on common models. The structure
#' used to create the data is similar as the \code{bamlss.formula}.
#'
#' @param formula List of the parameters, indicating how they should be computed.
#' similar to formula for \code{lm}, \code{glm}, \code{bamlss}, with the difference
#' that it included the coefficients and link function explicitly.
#'
#' @param generator Function to generate the response variables given the parameters
#'
#' @param n Number of observations to be simulated
#'
#' @param init_data Initial data including some variables to not been simulated.
#'
#' @param seed Seed to be defined with function \code{set.seed} to obtain reproducible
#' results
#'
#' @param extent Spatial extent for the simulation of coordinates when a spatial effect
#' is included.
#'
#' @return a \code{tibble} containing the simulated predictors, parameters and response
#' variable
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' f <- list(
#'   mean ~ I(5 + 0.5 * x1 + 0.1 * x2 + 0.7 * id1),
#'   sd ~ exp(x1)
#' )
#' (model.frame.sim(f, n = 10))
#' (data <- sim_model(f, n = 100))
#'
#' # Structure of the model
#' formula <- list(
#'   mean ~ I(age ^ 2) + fa(sex, beta = c(-1, 1)) + gp(s1, s2, beta = c(-1, 0)),
#'   sd ~ fa(sex, beta = c(1, 2))
#' )
#' idata <- data.frame(s1 = 1:10)
#' (datasim <- model.frame.sim(formula, idata = idata))
#' (datasim <- model.frame.sim(formula, n = 10))
#' (data <- model.frame.sim(formula, n = 10))
#' (model.response(data))
#' (datasim2 <- sim_model(formula, n = 10))
#' model.frame.sim(formula, n = 10) %>% model.response()
#'
#' @importFrom purrr map map_chr reduce
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#'
#' @export

sim_model <- function (formula = list(mean ~ I(1 + 2 * x1), sd ~ 1),
                       link = list(identity, exp), generator = rnorm,
                       n = nrow(init_data), q = 1, init_data = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  data <- model.frame.sim(formula, n = n, idata = init_data)
  data <- model.response(model_frame = data, link = link, generator = generator, q = 1)

  return(data)
}


#' @title Simulate response variable
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#'
#' f <- list(
#'   mean ~ I(5 + 0.5 * x1 + 0.1 * x2 + 0.7 * id1),
#'   sd ~ I(x1)
#' )
#' (model_frame <- model.frame.sim(f, n = 10))
#' (data <- model.response(model_frame, link = list(identity, exp)))
#'
#'
#' @export
model.response <- function (model_frame, formula = attr(model_frame, "formula"),
                            link = list(identity, exp), generator = rnorm,
                            q = 1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(model_frame)
  params <- purrr::map_chr(formula, ~ all.vars(.)[1])

  model_frame[params] <- purrr::map(formula, ~ .[[3]]) %>%
    purrr::map(~ eval(., model_frame)) %>%
    purrr::map(as.numeric) %>%
    purrr::map2(link, ~ .y(.x))

  model_frame["y"] <- do.call(generator, c(n = n, model_frame[params]))

  return(model_frame)
}

