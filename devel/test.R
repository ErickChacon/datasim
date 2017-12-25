
rm(list = ls())

# Structure of the model
formula <- list(
  mean ~ I(popo ^ 2) + mfe(vars1, beta = c(0.1, 0, 1)) + mgp(s1, s2, bla = exp(2)),
  # mean ~ mfe(vars1, beta = c(0.1, 0, 1)),
  sd ~ mfe(vars1, beta = c(-1, 1, 0.5))
)
generator = rnorm
n = 3000
# idata = NULL
seed = NULL
extent = 1
idata = data.frame(s1 = 1:n)
library(tidyverse)

datasim <- model.frame.sim(formula, idata = idata)
datasim
plot(datasim$s1, datasim$s2)

