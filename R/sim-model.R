
#' @title Simulate databases based in models
#'
#' @description
#' \code{sim_model} simulate a database based on common models. The structure
#' used to create the data is similar as the \code{bamlss.formula}.
#'
#' @param formula List of the parameters, indicating how they should be computed.
#' similar to formula for \code{lm}, \code{glm}, \code{bamlss}, with the difference
#' that it included the coefficients and link function explicitly.
#' @param link_inv A list of function representing the inverse link function for the
#' parameters.
#' @param generator Function to generate the response variables given the parameters
#' @param n Number of observations to be simulated
#' @param responses a numeric value indicating the number of response or a character
#' vector indicating the names of the response variables
#' @param init_data Initial data including some variables to not been simulated.
#' @param lm_syntax Optional logical argument to indicate if the usual syntax of
#' \code{formula.lm} should be used.
#' @param effects_save Optional logical argument to save or not generated random
#' effects
#' @param seed Seed to be defined with function \code{set.seed} to obtain reproducible
#' results
#'
#' @return A \code{tibble} containing the simulated predictors, parameters and response
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
#' (model_frame(f, n = 10))
#' (data <- sim_model(f, n = 100))
#'
#' # Structure of the model
#' formula <- list(
#'   mean ~ I(age ^ 2) + fa(sex, beta = c(-1, 1)),
#'   sd ~ fa(sex, beta = c(1, 2))
#' )
#' idata <- data.frame(s1 = 1:10)
#' (datasim <- model_frame(formula, idata = idata))
#' (datasim <- model_frame(formula, n = 10))
#' (data <- model_frame(formula, n = 10))
#' (model_response(data))
#' (datasim2 <- sim_model(formula, n = 10))
#' library(magrittr)
#' model_frame(formula, n = 10) %>% model_response()
#'
#' f <- list(
#'   mean ~ gp(list(s1), "exp_cor", list(phi = 0.05)),
#'   sd ~ I(1.6)
#' )
#' (data <- sim_model(f, n = 400))
#' plot(mean ~ s1, data)
#' plot(response ~ s1, data)
#'
#' formula <- list(
#' mean ~ I(0.5 * x1) : I(x2) + re(city, 1, 2),
#' sd ~ I(1)
#' )
#' data <- sim_model(formula, n = 10)
#'
#' @importFrom purrr map map_chr reduce
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#'
#' @export
sim_model <- function (formula = list(mean ~ I(0), sd ~ I(1)),
                       link_inv = replicate(length(formula), identity),
                       generator = rnorm, n = nrow(init_data),
                       responses = 1, init_data = NULL,
                       lm_syntax = TRUE,
                       effects_save = TRUE,  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  resp_list <- responses_check(responses)
  response_name <- resp_list$name
  response_labels <- resp_list$labels
  q <- resp_list$q

  data <- model_frame(formula, n = n, idata = init_data, q = q)

  if (lm_syntax == TRUE) {
    data <- model_response_lm(model_frame = data, link_inv = link_inv,
                               generator = generator, responses = responses,
                               effects_save = effects_save)
  } else {
    data <- model_response(model_frame = data, link_inv = link_inv,
                           generator = generator, responses = responses)

  }

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
#' @param formula List of the parameters, indicating how they should be computed.
#' similar to formula for \code{lm}, \code{glm}, \code{bamlss}, with the difference
#' that it included the coefficients and link function explicitly.
#' @param n Number of observations to be simulated
#' @param idata Initial data including some variables to not been simulated.
#' @param q number of response variables
#' @param seed Seed to be defined with function \code{set.seed} to obtain reproducible
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' # Structure of the model
#' formula <- list(
#'   mean ~ I(age ^ 2) + fa(sex, beta = c(-1, 1)),
#'   sd ~ fa(sex, beta = c(1, 2))
#' )
#' idata <- data.frame(s1 = 1:10)
#' (datasim <- model_frame(formula, idata = idata))
#' (datasim <- model_frame(formula, n = 10))
#'
#' formula <- list(
#'   mean ~ I(5) + I(0.5 * x1) + I(0.1 * x2),
#'   sd ~ I(2) + I(x)
#' )
#'
#' # Structure of the model
#' formula <- list(
#'   mean ~ fa(sex, beta = c(-1,1)),
#'   sd ~ I(0)
#' )
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map reduce map_chr
#' @importFrom dplyr select mutate arrange group_by mutate ungroup bind_cols n
#' @importFrom tidyr unnest
#' @importFrom stats setNames
#'
#' @export
model_frame <- function (formula, n = nrow(idata) / q, idata = NULL,
                         q = 1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Effects that can generate covariates
  generators <- c("mfa", "fa", "mre", "re", "mgp", "gp")

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
  generators_order <- rev(unique(c(generators, effects$type)))
  covariates <- tidyr::unnest(dplyr::select(effects, -call)) %>%
    dplyr::mutate(type_order = match(type, generators_order)) %>%
    dplyr::arrange(type_order) %>%
    dplyr::group_by(covs) %>%
    dplyr::mutate(
      rep = 1:dplyr::n(),
      generate = "none"
      ) %>%
    dplyr::ungroup()

  # Establish how and which covariates to generate
  if (is.null(idata)) {
    idata <- tibble::tibble(id = rep(1:n, q))
  } else {
    if (nrow(idata) == n & q > 1) {
      idata <- do.call("rbind", replicate(q, idata, simplify = FALSE))
    }
    idata$id <- rep(1:n, q)
  }
  covariates <- within(covariates, {
    generate[!(type %in% generators)] <- "gaussian"
    generate[type %in% generators] <- "generator"
    generate[covs %in% names(idata)] <- "none"
    generate[rep != 1] <- "none"
  })

  # Which call effects should be executed and initialize data
  gener_exec <- unique(with(covariates, id[generate == "generator"]))
  data <- tibble::tibble(id = rep(1:n, q))

  for (i in gener_exec) {

    # Replace variables that need to be simulated to NA in the original call
    aux_covs <- covariates %>%
      dplyr::filter(rep == 1, covs %in% effects$covs[[i]])
    aux_covs_na <- purrr::map(aux_covs$generate, ~ ifelse(. == "generator", NA, .)) %>%
      stats::setNames(aux_covs$covs)
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
  data[gauss_names] <- replicate(length(gauss_names), rep(rnorm(n), q)) %>%
    as.data.frame()

  # Covariates order
  covariates_order <- c("id", unique(purrr::reduce(effects$covs, c)))
  data <- data[intersect(covariates_order, names(data))]

  # Join simulated data with provided data
  data <- dplyr::select(data, -id)
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
#' @param model_frame A \code{tibble} containing all the covariates to be used.
#' Usually this is obtained as the output of the \code{model_frame} function.
#' @param formula An optional list of formulas to simulate the parameters of the
#' response variables.
#' @param link_inv A list of function representing the inverse link function for the
#' parameters.
#' @param generator Function to generate the response variables given the parameters
#' @param responses character vector indicating the names of the response variables
#' @param seed Seed to be defined with function \code{set.seed} to obtain reproducible
#' results
#'
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
#' (model_fr <- model_frame(f, n = 10))
#' (data <- model_response(model_fr, link = list(identity, exp)))
#'
#' beta0 <- c(-1, 1)
#' # Structure of the model
#' formula <- list(
#'   mean ~ fa(sex, beta = get("beta0")),
#'   sd ~ I(0)
#' )
#' (model_fr <- model_frame(formula, n = 10))
#' (data <- model_response(model_fr, link = list(identity, exp)))
#'
#' formula <- list(
#'   mean ~ mfe(x, beta = 1:2),
#'   sd ~ mfe(x1, beta = 1:2)
#' )
#' (model_fr <- model_frame(formula, n = 10))
#' #(data <- model_response(model_fr, link_inv = list(identity, exp), responses = 1:2))
#' #(data <- sim_model(formula, link_inv = list(identity, exp), n = 10, responses = 1:2))
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map_chr map2
#' @importFrom dplyr left_join
#'
#' @export
model_response <- function (model_frame, formula = attr(model_frame, "formula"),
                            link_inv = replicate(length(formula), identity),
                            generator = rnorm,
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
  model_frame <- dplyr::left_join(model_frame, tibble::as_tibble(params_ls), by = "id")

  return(model_frame)
}


#' @title Simulate response variable (testing function)
#'
#' @description
#' \code{function} description.
#'
#' @param model_frame A \code{tibble} containing all the covariates to be used.
#' Usually this is obtained as the output of the \code{model_frame} function.
#' @param formula An optional list of formulas to simulate the parameters of the
#' response variables.
#' @param link_inv A list of function representing the inverse link function for the
#' parameters.
#' @param generator Function to generate the response variables given the parameters
#' @param responses a numeric value indicating the number of response or a character
#' vector indicating the names of the response variables
#' @param effects_save Optional logical argument to save or not generated random
#' effects
#' @param seed Seed to be defined with function \code{set.seed} to obtain reproducible
#' results
#'
#' @examples
#'
#' formula <- list(
#' mean ~ I(0.5 * x1) : I(x2) + re(city, 1, 2),
#' sd ~ I(1)
#' )
#' model_frame <- model_frame(formula, n = 10)
#' model_response_lm(model_frame)
#'
#'
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map_chr map2
#' @importFrom dplyr left_join
#'
#' @export

model_response_lm <- function (model_frame, formula = attr(model_frame, "formula"),
                            link_inv = replicate(length(formula), identity),
                            generator = rnorm,
                            responses = c("response"),
                            effects_save = FALSE,
                            seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Get dimensions and parameters
  nq <- nrow(model_frame)
  n <- length(unique(model_frame$id))
  q <- nq / n
  params <- purrr::map_chr(formula, ~ all.vars(.)[1])
  nterms <- purrr::map(formula, ~ length(attr(terms(.), "variables")) - 2)

  # Get effects details from formula
  effects <- tibble::tibble(
    params = rep(params, nterms),
    call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
      purrr::map(~ .[c(-1, -2)]) %>%
      purrr::reduce(c),
    # call_str = format(call),
    call_str = purrr::map_chr(call, ~ deparse(., width.cutoff = 300)),
    # call_str = purrr::map(formula, ~ labels(terms(.))) %>%
    #   purrr::reduce(c),
    id = 1:length(call),
    covs = purrr::map(call, all.vars),
    type = purrr::map_chr(call, ~ as.character(.x[[1]])),
    ncovs = purrr::map_int(covs, length)
    )

  # Remove intercepts from formula
  param_form <- effects %>%
    dplyr::filter(ncovs == 0) %>%
    dplyr::left_join(tibble::tibble(params), ., by = "params") %>%
    dplyr::mutate(
      offset = purrr::map_dbl(call, ~ ifelse(is.null(.), 0, eval(.))),
      update_formula = purrr::map(call, ~ substitute( ~ . - x, list(x = .))),
      formula_new = purrr::map2(formula, update_formula, ~ update(.x[-2], .y)),
      x = purrr::map(formula_new, ~ model.frame(., model_frame)),
      x_mat = purrr::map2(formula_new, x, ~ model.matrix(.x, .y)),
      value = purrr::map2(x_mat, offset,
        ~ if(ncol(.x) > 1) rowSums(.x) - 1 + .y else .y)
      )

  # Select effects to be exported
  generators2save <- c("mre", "re", "mgp", "gp")
  param_data <- param_form %>%
    dplyr::mutate(
      eff_show = purrr::map(params,
        ~ subset(effects, params == . & type %in% generators2save,
                 select = call_str, drop = TRUE)),
      x_show = purrr::map2(x, eff_show, ~ dplyr::select(.x, .y)),
      x_show = purrr::map2(x_show, params,
        purrr::possibly(~ setNames(.x, paste0(names(.x), ".", .y)), NULL))
      )
  param_data <- dplyr::bind_cols(param_data$x_show)
  names(param_data) <-
    gsub("^(\\w+)\\((\\w+).*?(\\w+)$", "\\1.\\2.\\3", names(param_data)) %>%
    make.names(unique = TRUE)
  # if (ncol(param_data) > 0) {
  #   param_data$id <- rep(1:n,
  #   names(param_data) <- make.names(names(param_data), unique = TRUE) %>%
  #     gsub("\\.+", "\\.", .)
  # } else {
  #   param_data <- tibble::tibble(id = 1:10)
  # }


  # Compute parameters in a new list
  params_ls <- list()
  params_ls[params] <- param_form$value %>%
    purrr::map2(link_inv, ~ .y(.x))

  # Simulate response variables on the list
  resp_list <- responses_check(responses)
  response_name <- resp_list$name
  response_labels <- resp_list$labels
  params_ls[response_name] <- list(do.call(generator, c(n = n * q, params_ls[params])))
  response_label <- paste0(response_name, "_label")
  if (q > 1) params_ls[response_label] <- list(rep(response_labels, each = n))
  # params_ls$id <- rep(1:n, q)

  # Join covariates, with parameters and response variables
  if (effects_save == TRUE) {
    # model_frame <- dplyr::left_join(model_frame, param_data, by = "id")
    model_frame <- dplyr::bind_cols(model_frame, param_data)
  }
  model_frame <- dplyr::bind_cols(model_frame, tibble::as_tibble(params_ls))

  return(model_frame)
}

responses_check <- function (responses = "response") {
  if (length(responses) == 1) {
    if (is.character(responses)) {
      q <- 1
      name <- responses
      labels <- NULL
    } else {
      q <- responses
      name <- "response"
      labels <- 1:q
    }
  } else if (length(responses) > 1) {
    q <- length(responses)
    name <- "response"
    if (is.character(responses)) {
      labels <- responses
    } else {
      labels <- 1:q
    }
  }
  return(list(q = q, name = name, labels = labels))
}
