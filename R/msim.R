
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



beta0 <- -1:1
data <- tibble(id = 1:10)
(data[c("s1", "s2")] <- gp(s1 = NA, s2 = NA, beta = -1:1, size = 10))
eval(quote(gp(s1, s2, beta = beta0)), data)


(s1 <- x[[1]])
(s2 <- x[[2]])
gp(s1, s2, beta = -1:0)
