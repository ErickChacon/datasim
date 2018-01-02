# library(datasim)
#
# f <- list(
#   prob ~ mfa(ones, beta = rbind(1:3)) +
#     re(id, sigma = 1, groups = 6, q = 3):mfa(ones, beta = rbind(2:4)),
#   size ~ I(1)
#   )
# formula <- f
# idata <- NULL
# n <- 6
# q <- 3
#
# data_model <- sim_model(formula = f,
#                         # link_inv = list(psych::logistic, identity),
#                         link_inv = list(psych::logistic, identity),
#                         generator = rbinom,
#                         responses = 3,
#                         n = 6)
# data_model
#
# # formula <- list(
# #   mean ~ I(0.5 * x1) : I(x2) + re(city, 1, 2),
# #   sd ~ I(1)
# # )
#
# # model_frame <- model_frame(f, n = 5, q = 3)
# data <- model_frame(f, n = 5, q = 3)
# formula = attr(model_frame, "formula")
# link_inv = list(psych::logistic, identity)
# generator = rbinom
# responses = 3
# effects_save = FALSE
#
# data <- model_response_lm(model_frame = data, link_inv = link_inv,
#                            generator = generator, responses = responses,
#                            effects_save = effects_save)
#
# # Get effects details from formula
# effects <- tibble::tibble(
#   call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
#     purrr::map(~ .[c(-1, -2)]) %>%
#     purrr::reduce(c),
#   id = 1:length(call),
#   covs = purrr::map(call, all.vars),
#   type = purrr::map_chr(call, ~ as.character(.x[[1]]))
# )
#
# param_formula <- tibble::tibble(
#   formula = formula,
#   params = purrr::map_chr(formula, ~ all.vars(.)[1]),
#   call = purrr::map(formula, ~ as.list(attr(terms(.), "variables"))) %>%
#     purrr::map(~ .[c(-1, -2)]),
#   call_str = purrr::map(call, format),
#   covs = purrr::map(call, all.vars),
#   # %>%
#   #     purrr::reduce(c),
#   )
# param_formula #%>% tidyr::unnest(call_str)
# param_formula$call[[1]]
#
#   # Compute parameters in a new list
#   params_ls <- list()
#   params_ls[params] <- purrr::map(formula, ~ .[[3]]) %>%
#     purrr::map(~ eval(., model_frame)) %>%
#     purrr::map(as.numeric) %>%
#     purrr::map2(link_inv, ~ .y(.x))
#
#
