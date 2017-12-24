
rm(list = ls())

model.frame.sim <- function (formula, n = nrow(idata), idata = NULL, seed = NULL) {

  generators <- c("mfe", "mgp")
  effects <- tibble(
    call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
      purrr::map(~ .[c(-1, -2)]) %>%
      purrr::reduce(c),
    id = 1:length(call),
    covs = purrr::map(call, all.vars),
    type = purrr::map_chr(call, ~ as.character(.x[[1]]))
  )
  effects
  covariates <- unnest(dplyr::select(effects, -call)) %>%
    group_by(covs, type) %>%
    mutate(
      rep = 1:n(),
      generate = "none"
      ) %>%
    ungroup()
  covariates <- within(covariates, {
    generate[!(type %in% generators)] <- "gaussian"
    generate[type %in% generators & rep == 1] <- "generator"
    generate[covs %in% names(idata)] <- "none"
  })

  covariates

  gener_exec <- unique(with(covariates, id[generate == "generator"]))

  data <- tibble(id = 1:n)

  for (i in gener_exec) {

    # Replace variables that need to be simulated to NULL in the original call
    aux_covs <- covariates %>%
      filter(rep == 1, covs %in% effects$covs[[i]])
    aux_covs_na <- map(aux_covs$generate, ~ ifelse(. == "generator", NA, .)) %>%
      setNames(aux_covs$covs)
    aux_call <- do.call('substitute', list(effects$call[[i]], aux_covs_na))
    aux_call

    # Adding sample size to call
    aux_call <- as.list(aux_call)
    aux_call[["size"]] <- n
    aux_call <- as.call(aux_call)
    aux_call

    # Generate required covariates
    aux_names <- filter(aux_covs, generate == "generator") %>% pull("covs")
    data[aux_names] <- eval(aux_call)
    data

  }

  data <- dplyr::bind_cols(idata, dplyr::select(data, -id))
  return(data)
}


mfe <- function (x = NA, beta, labels = 1:length(beta), size = NULL) {
  if (is.na(x)) {
    output <- factor(sample(labels, size = size, replace = TRUE), labels = labels)
  } else {
    names(beta) <- levels(x)
    output <- as.numeric(beta[x])
  }
  return(output)
}

mgp <- function (..., beta, size = NULL) {
  vars <- list(...)
  nvar <- purrr::map(vars, is.na) %>% do.call(sum, .)
  if (nvar > 0) {
    output <- replicate(nvar, list(runif(size)))
  } else {
    output <- s1 * beta[1] + s2 * beta[2]
  }
  return(output)
}


# Structure of the model
formula <- list(
  mean ~ I(popo ^ 2) + mfe(vars1, beta = c(0.1, 0, 1)) + mgp(s1, s2, bla = exp(2)),
  # mean ~ mfe(vars1, beta = c(0.1, 0, 1)),
  sd ~ mfe(vars1, beta = c(-1, 1, 0.5))
)
generator = rnorm
n = 10
# idata = NULL
seed = NULL
extent = 1
idata = data.frame(s1 = 1:10)
library(tidyverse)

datasim <- model.frame.sim(formula, idata = idata)
datasim

