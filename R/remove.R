# library(datasim)
# formula <- list(
#   mean ~ I(0.5 * x1) : I(x2) + re(city, 1, 2),
#   sd ~ I(1)
# )
# model_frame <- model_frame(formula, n = 10)
# formula = attr(model_frame, "formula")
# link_inv = list(identity, exp)
# generator = rnorm
# responses = c("response")
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
