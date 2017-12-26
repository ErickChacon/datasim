
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
                       link_inv = list(identity, exp), generator = rnorm,
                       n = nrow(init_data), responses = c("response"), init_data = NULL,
                       seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  data <- model.frame.sim(formula, n = n, idata = init_data)
  data <- model.response(model_frame = data, link_inv = link_inv, generator = generator,
                         responses = responses)

  return(data)
}


#' @title Simulate covariates based on model formula
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
#' # Structure of the model
#' formula <- list(
#'   mean ~ I(age ^ 2) + fa(sex, beta = c(-1, 1)) + gp(s1, s2, beta = c(-1, 0)),
#'   sd ~ fa(sex, beta = c(1, -1))
#' )
#' idata <- data.frame(s1 = 1:10)
#' (datasim <- model.frame.sim(formula, idata = idata))
#' (datasim <- model.frame.sim(formula, n = 10))
#'
#' formula <- list(
#'   mean ~ I(5 + 0.5 * x1 + 0.1 * x2 + 0.7 * id),
#'   sd ~ I(x1)
#' )
#'
#' beta0 <- c(-1, 1)
#' # Structure of the model
#' formula <- list(
#'   mean ~ fa(sex, beta = beta0),
#'   sd ~ I(0)
#' )
#'
#' @export
model.frame.sim <- function (formula, n = nrow(idata), idata = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Effects that can generate covariates
  generators <- c("fa", "gp")

  # Get effects details from formula
  effects <- tibble::tibble(
    call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
      purrr::map(~ .[c(-1, -2)]) %>%
      purrr::reduce(c),
    id = 1:length(call),
    covs = purrr::map(call, all.vars),
    type = purrr::map_chr(call, ~ as.character(.x[[1]]))
  )

  # Get covariates details from effects
  generators_order <- rev(unique(c(generators_order, effects$type)))
  covariates <- unnest(dplyr::select(effects, -call)) %>%
    mutate(type_order = match(type, generators_order)) %>%
    arrange(type_order) %>%
    group_by(covs) %>%
    mutate(
      rep = 1:n(),
      generate = "none"
      ) %>%
    ungroup()

  # Establish how and which covariates to generate
  covariates <- within(covariates, {
    generate[!(type %in% generators)] <- "gaussian"
    generate[type %in% generators] <- "generator"
    generate[covs %in% names(idata)] <- "none"
    generate[rep != 1] <- "none"
  })

  # Which call effects should be executed and initialize data
  gener_exec <- unique(with(covariates, id[generate == "generator"]))
  data <- tibble::tibble(id = 1:n)

  for (i in gener_exec) {

    # Replace variables that need to be simulated to NA in the original call
    aux_covs <- covariates %>%
      dplyr::filter(rep == 1, covs %in% effects$covs[[i]])
    aux_covs_na <- purrr::map(aux_covs$generate, ~ ifelse(. == "generator", NA, .)) %>%
      setNames(aux_covs$covs)
    aux_call <- do.call('substitute', list(effects$call[[i]], aux_covs_na))

    # Adding sample size to call
    aux_call <- as.list(aux_call)
    aux_call[["size"]] <- n
    aux_call <- as.call(aux_call)

    # Generate required covariates
    aux_names <- dplyr::filter(aux_covs, generate == "generator") %>%
      dplyr::pull("covs")
    data[aux_names] <- eval(aux_call)

  }

  # Generate covariates with Gaussian distribution
  gauss_names <- with(covariates, covs[generate == "gaussian"])
  data[gauss_names] <- as.data.frame(replicate(length(gauss_names), rnorm(n)))

  # Covariates order
  covariates_order <- c("id", unique(purrr::reduce(effects$covs, c)))
  data <- data[intersect(covariates_order, names(data))]

  # Join simulated data with provided data
  data <- dplyr::bind_cols(idata, data)
  attr(data, "formula") <- formula

  return(tibble::as_tibble(data))
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
#' beta0 <- c(-1, 1)
#' # Structure of the model
#' formula <- list(
#'   mean ~ fa(sex, beta = beta0),
#'   sd ~ I(0)
#' )
#' (model_frame <- model.frame.sim(formula, n = 10))
#' (data <- model.response(model_frame, link = list(identity, exp)))
#'
#' formula <- list(
#'   mean ~ mfe(x, beta = 1:2),
#'   sd ~ mfe(x1, beta = 0:1)
#' )
#' (model_frame <- model.frame.sim(formula, n = 10))
#' (data <- model.response(model_frame, responses = 1:2))
#' (data <- sim_model(formula, n = 10, response = 1:2))
#'
#' @export
model.response <- function (model_frame, formula = attr(model_frame, "formula"),
                            link_inv = list(identity, exp), generator = rnorm,
                            responses = c("response"), seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Get dimensions and parameters
  n <- nrow(model_frame)
  q <- length(responses)
  params <- purrr::map_chr(formula, ~ all.vars(.)[1])

  # Compute parameters in a new list
  params_ls <- list()
  params_ls[params] <- purrr::map(formula, ~ .[[3]]) %>%
    purrr::map(~ eval(., model_frame)) %>%
    purrr::map(as.numeric) %>%
    purrr::map2(link_inv, ~ .y(.x))

  # Simulate response variables on the list
  response <- ifelse(q > 1, "response", responses)
  params_ls[response] <- list(do.call(generator, c(n = n * q, params_ls[params])))
  if (q > 1) params_ls$response_label <- rep(responses, each = n)
  params_ls$id <- rep(1:n, q)

  # Join covariates, with parameters and response variables
  model_frame <- dplyr::left_join(model_frame, as_tibble(params_ls), by = "id")

  return(model_frame)
}

