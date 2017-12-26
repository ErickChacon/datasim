
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
#s' @export
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
  covariates_order <- unique(purrr::reduce(effects$covs, c))
  data <- data[intersect(covariates_order, names(data))]

  # Join simulated data with provided data
  data <- dplyr::bind_cols(idata, data)
  attr(data, "formula") <- formula

  return(tibble::as_tibble(data))
}

