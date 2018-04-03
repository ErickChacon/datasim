
library(datasim)
n <- 100
q <- 3

    # mfe(x, beta = c(0.1, 0, 1)),
f <- list(
  # mean ~ mfa(x, beta = cbind(1:2, 2:3, 0:1)),
  # mean ~ mfe(x, beta = c(0.1, 0, 1)),
  mean ~ mfa(x, beta = cbind(1:2, 2:3, 0:1)) + mfe(y, beta = c(0.1, 0, 1)),
  sd ~ I(1)
  )

  x <- mfa(beta = cbind(1:2, 2:3, 0:1))
x <- rep(1, 15)
  mfa(x, beta = cbind(1:2, 2:3, 0:1))

  mfe(x, beta = c(0.1, 0, 1))

(model_frame <- model_frame(f, n = 5, q = 3))

  (x <- mfa(beta = cbind(1:2, 2:3, 0:1), size = 5))
  mfa(x, beta = cbind(1:2, 2:3, 0:1))

  (mfe(x <- rep(1:5,3), beta = c(0.1, 0, 1)))

data_model <- sim_model(formula = f, n = 100, seed = 1, responses = 3)

knitr::kable(head(data_model, 300))

roxygen2::roxygenize(getwd(), roclets = c("collate", "namespace", "rd"))


# # model_frame <- model_frame(f, n = 5, q = 3)
# model_frame <- model_frame(f, n = n, q = q)
# formula = attr(model_frame, "formula")
#
# link_inv = list(identity, identity)
# generator = rbinom
# responses = q
# effects_save = TRUE
#
#
# data_long <- sim_model(formula = f,
#                         link_inv = list(identity, identity),
#                         generator = rnorm,
#                         responses = c("a", "b", "c"),
#                         n = n,
#                         seed = 2)
#
#
#
