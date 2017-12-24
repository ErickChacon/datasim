
# Covariance parameters
n <- 100
q <- 2
var <- sqrt(diag(c(4, 4)))
A <- matrix(c(1, - 0.8, 0, 0.6), nrow = 2)
variance <- var %*% tcrossprod(A) %*% var
nugget <- diag(0, q)
phi <- rep(1 / 0.08, q)

# Structure of the model
formula <- list(
  mean ~ psych::logistic(
    mgp(s1, s2, "exponential", get("variance"), get("nugget"), get("phi"))),
  sd ~ 1
)

# Simulate data based on formula
library(tidyr)
library(dplyr)
data <- msim_model(formula, generator = rnorm, n = n, extent = 2, seed = 1)
data_long <- gather(data, yname, yval, matches("^y[0-9]+"))

# Plot the observed realization
library(ggplot2)
spgg <- ggplot(data_long, aes(s1, s2, size = yval, col = yval)) +
  geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  facet_wrap(~ yname)
print(spgg)




mfe(x = rnorm(10), beta = c(0.1, 0, 1))

# -----

mfe <- function (x, beta) {
  as.numeric(matrix(x) %*% matrix(beta, nrow = 1))
}

# Structure of the model
formula <- list(
  mean ~ mfe(vars1, beta = c(0.1, 0, 1)) + mgp(s1, s2, bla = exp(2)),
  # mean ~ mfe(vars1, beta = c(0.1, 0, 1)),
  sd ~ 1
)
generator = rnorm
n = 100
init_data = NULL
seed = NULL
extent = 1
init_data = data.frame(s1 = 1:10)
library(tidyverse)

msim_model <- function (formula, generator = rnorm, n = 100, init_data = NULL,
                        seed = NULL, extent = 1) {

  if (!is.null(seed)) set.seed(seed)

  formula <- mean ~ mfe(vars1, beta = c(0.1, 0, 1)) + mgp(s1, s2, bla = exp(2)) + re(pla)
  formula1 <- mean ~ mfe(vars1, beta = c(0.1, 0, 1))
  as.character(formula)
  as.character(formula1)
  as.character(formula[[3]])
  as.character(formula1[[3]][[3]])
  class(formula1[[3]][[3]])
  all.vars(formula1[[3]])
  class(formula1[[3]])

  plopa <- formula[[1]][[3]]
  class(formula[[1]][[3]])
  all.vars(formula[[1]], functions = TRUE)
  all.vars(formula[[1]])
  # all.names(formula[[1]])

  formula[[3]]
  termsss <- terms(formula)
  termss <- labels(terms(formula1))
  term0 <- parse(text = termss[1])
  all.vars(term0)

  all.vars(parse(text = termss[1]))
  all.names(parse(text = termss[1]))
  all.vars(parse(text = termss[2]))
  all.vars(parse(text = termss[3]))

  as.character(formula[[3]])
  formula[[3]][[2]]
  formula[[3]][[3]]

  formula1[[3]]
  as.character(formula[[1]])
  plopa[[1]]
  plopa[[3]]

  as.character(formula[[1]][[3]])
  form_text <- format(formula[[1]][[3]])
  form_text

  aux <- regmatches(form_text, gregexpr("mfe.*?,", form_text))[[1]]
  sub("mfe\\(([[:alnum:]]+),", "\\1", aux)

  aux <- regmatches(form_text, gregexpr("mgp.*?\\)", form_text))[[1]]
  sub("mfe\\(([[:alnum:]]+),", "\\1", aux)

  params <- purrr::map_chr(formula, ~ all.vars(.)[1])
  predictors <- purrr::map(formula, ~ all.vars(.)[-1]) %>%
    purrr::reduce(c) %>% unique()
  p <- length(predictors)

  # Identify predictors that need to be simulated
  init_pred <- names(init_data)
  pred2sim <- setdiff(predictors, init_pred)
  pred_sp <- grep("^s[0-9]+$", pred2sim, value = TRUE)
  pred2sim <- setdiff(pred2sim, pred_sp)
  p2 <- length(pred2sim)
  p_sp <- length(pred_sp)

  init_pred
  pred2sim
  pred_sp

  predictors_data <- data.frame(predictors = predictors)
  predictors_data$init <- predictors_data$predictors %in% init_pred
  predictors_data$sp <- predictors_data$predictors %in% pred_sp
  predictors_data$gaussian <- apply(subset(predictors_data, select = init:sp), 1, sum) == 0
  predictors_data

  

  # Simulate only required predictors
  if (p2 > 0) {
    data <- matrix(rnorm(n * p2), nrow = n) %>%
      tibble::as_tibble() %>%
      setNames(pred2sim)
    init_data <- dplyr::bind_cols(init_data, data)
  }
  if (p_sp > 0) {
    data <- matrix(runif(n * p_sp) * extent, nrow = n) %>%
      tibble::as_tibble() %>%
      setNames(pred_sp)
    init_data <- dplyr::bind_cols(init_data, data)
  }

  init_data$id <- 1:nrow(init_data)

  # Evaluate parameters in a tibble
  params_ls <- list()
  params_ls[params] <-
    
    purrr::map(formula, ~ .[[3]]) %>%
    purrr::map(~ eval(., init_data))
  params_ls <- tibble::as_tibble(params_ls)

  # Obtain dimensions
  nq <- max(purrr::map_int(params_ls, length))
  q <- round(nq / n)

  # Simulate multivariate process and reshape
  params_ls$y <- do.call(generator, c(n = nq, params_ls[params]))
  params_ls$number <- rep(1:q, each = n)
  params_ls$id <- rep(1:n, q)
  # varname <- value <- number <- id <- NULL
  params_ls <- params_ls %>%
    tidyr::gather("varname", "value", - number, - id) %>%
    dplyr::mutate(varname = paste0(varname, number)) %>%
    dplyr::select(- number) %>%
    tidyr::spread(varname, value)

  # Organize and joint final simulated dataset
  init_data <- dplyr::left_join(init_data, params_ls, by = "id")

  return(init_data)
}


# Structure of the model
formula <- list(
  mean ~ I(popo ^ 2) + mfe(vars1, beta = c(0.1, 0, 1)) + mgp(s1, s2, bla = exp(2)),
  # mean ~ mfe(vars1, beta = c(0.1, 0, 1)),
  sd ~ mfe(vars1, beta = c(-1, 1, 0.5))
)
generator = rnorm
n = 10
init_data = NULL
seed = NULL
extent = 1
init_data = data.frame(s1 = 1:10)
library(tidyverse)

model.frame.sim <- function (formula, n = 100, init_data = NULL, seed = NULL) {

  formula

  params <- purrr::map_chr(formula, ~ all.vars(.)[1])
  predictors <- purrr::map(formula, ~ all.vars(.)[-1]) %>%
    purrr::reduce(c) %>% unique()
  p <- length(predictors)

  generators <- c("mfe", "mgp")
  effects <- tibble(
    call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
      purrr::map(~ .[c(-1, -2)]) %>%
      purrr::reduce(c),
    id = 1:length(call),
    covs = purrr::map(call, all.vars),
    type = purrr::map_chr(call, ~ as.character(.x[[1]])),
  )
  effects
  covariates <- unnest(dplyr::select(effects, -call)) %>%
    group_by(covs, type) %>%
    mutate(
      rep = 1:n(),
      generate = NA
      )
  covariates <- within(covariates,
    {
      generate[!(type %in% generators)] <- "gaussian"
      generate[type %in% generators & rep == 1] <- "generator"
      generate[covs %in% names(init_data)] <- "none"
    }
  )
  covariates



  formula[[1]][[3]]
  purrr::map(formula, ~ labels(terms(.))) %>%
    purrr::reduce(c) %>%
    purrr::map(~ parse(text = .)) %>%
    purrr::map(~ as.list(quote, .))

  text = "a + b"
  do.call(quote, list(text))

  %>% purrr::map(as.call)

  plop <- as.list(attr(terms(formula[[2]]), "variables"))
  plop[c(-1, -2)]
  plop[[1:2]]

  # plop[[4]]
  plop[[1]] <- plop[[2]] <- NULL
  plop

  class(as.list(attr(terms(formula[[1]]), "variables"))[[3]])
  str(attr(terms(formula[[1]]), "factors"))
  rownames(attr(terms(formula[[1]]), "factors"))
  all.vars(formula[[1]])
  all.names(formula[[1]])

  termsss <- terms(formula[[1]])
  termss <- labels(terms(formula[[1]]))
  term0 <- parse(text = termss[1])
  all.vars(term0)


}

mfe.covs <- function () {

}

x <- factor(rep(c(1, 3, 7), each = 5))
y <- 1:15
model.frame(~ x + y)
beta <- c(0.1, 0, -0.3)

mfe <- function (x = NA, beta, labels = 1:length(beta), size = NULL) {
  if (is.na(x)) {
    output <- factor(sample(labels, size = size, replace = TRUE), labels = labels)
  } else {
    names(beta) <- levels(x)
    output <- as.numeric(beta[x])
  }
  return(output)
}
x <- mfe(NA, beta = c(-1, 1), size = 10, labels = c("male", "female"))
x

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

s <- mgp(s1 = NA, beta = c(-1, 1), size = 10)
s

 mgp(s1 = NA, 1, beta = c(-1, 1), size = 10)

x2 <- list(x, x)
data <- data.frame(id = 1:10)
data[c("x1", "x2")] <- x2
data["plopijl"] <- list(x)

# update call
blas <- quote(mfe(x, beta = c(-1, 1)))
blac <- as.list(blas)
blac$size = 10
blac <- as.call(blac)
eval(blac)

call(blac)

update(blas)
match.call(mfe, quote(mfe(x, beta = c(-1, 1))))

match.call(mfe, call("mfe", NULL, c(-1,1)), expand.dots = TRUE)

test <- call("mfe", NULL, c(-1,1))
(as.character(test))

ca <- substitute(sum(x,y))
class(ca)
match.call(sum, ca)

expr <- substitute(plop <- mfe(x, beta = c(-1, 1), labels = c("male", "female")))
match.call(get, call("get", "abc", i = FALSE, p = 3))

expr <- 'mfe(x, beta = c(-1, 1), labels = c("male", "female"))'
expr <- parse(text = sub("\\)$", ", size = 100)", expr))
bla <- as.call(expr)
x <- eval(do.call('substitute', list(expr[[1]], list(x = NULL))))
x
eval(expr)

expr <- expression(mfe(x, beta = c(-1, 1), size = 10))
x <- eval(do.call('substitute', list(expr[[1]], list(x = NULL))))
x
eval(expr)


# quote(cos(x))
# substitute(expr, list(x = NULL))




