#
# library(datasim)
#
# rm(list = ls())
# formula <- list(
#   mean ~ I(5) +
#     gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02), sigma2 = 1),
#   sd ~ I(1)
# )
# n <- 10
# # idata <- NULL
# # q <- 1
# # seed <- NULL
# model_frame <- model_frame(formula, n = n)
# link_inv <- replicate(length(formula), identity)
# generator <- rnorm
# responses = c("response")
# effects_save = TRUE
# seed = NULL
#
# data <- sim_model(formula = f,
#                   link_inv = list(pnorm, identity),
#                   generator = rbinom,
#                   n = n,
#                   seed = 1)
#
#
# mygp <- gp(list(data1$s1, data1$s2), cor.model = "exp_cor",
#            cor.params = list(phi = 0.02), sigma2 = 5)
# var(mygp)
#
# # data_model <- sim_model(formula = f, n = 500, seed = 1)
#
# (x <- gp(list(s1 = NA), cor.model = "exp_cor", cor.params = list(phi = 0.05), sigma2 = 2,
#          size = 10))
# (s1 <- x[[1]])
# # Simulate and plot the realization of a Gaussian process
# (y <- gp(list(s1), cor.model = "exp_cor", cor.params = list(phi = 0.05), sigma2 = 20))
# var(y)
#
# rm(list = ls())
# formula <- list(
#   mean ~ I(0.5 * x1) : I(x2) + re(city, 1, 2),
#   sd ~ I(1)
# )
# model_frame <- model_frame(formula, n = 10)
# # model_response_lm(model_frame)
# link_inv <- replicate(length(formula), identity)
# generator <- rnorm
# responses = c("response")
# effects_save = TRUE
# seed = NULL
#
#
# form <- y ~ gp(list(s1, s2),
#                cor.model = "exp_cor",
#                cor.params = list(phi = 0.02),
#                sigma2 = 1)
#
# # terms(form)
# # deparse(form[[3]], width.cutoff = 300)
# # format(form[[3]], trim = FALSE)
#
#
