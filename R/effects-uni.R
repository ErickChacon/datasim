
#' @title Univariate factor effects
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param x factor vector to evaluate the effect
#' @param beta numeric vector of the effects for each level of \code{x}
#' @param labels character vector to name the \code{levels} of \code{x}
#' @param size numeric value to simulate the covariate \code{x}
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' (x <- fa(beta = -1:1, size = 10))
#' fa(x, beta = -1:1)
#' fa(x, beta = 1:3)
#'
#' @export
fa <- function (x, beta, labels = 1:length(beta), size = NULL) {
  if (!is.null(size)) {
    output <- factor(sample(labels, size = size, replace = TRUE), levels = labels)
  } else {
    names(beta) <- levels(x)
    output <- as.numeric(beta[x])
  }
  return(output)
}

#' @title Univariate random effects
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param x factor vector to evaluate the random effect
#' @param sigma variance of the random effect
#' @param groups character vector to name the \code{levels} of \code{x} or a numeric
#' value indicating the number of groups of \code{x}
#' @param size numeric value to simulate the covariate \code{x}
#' @param q number of response varables to replicate
#' @param replace An optional argument to simulate an independent random effect with
#' one repetition (e.g. \code{mre(groups = 100, size = 100, replace = FALSE)})
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' (x <- re(groups = 5, size = 15))
#' re(x, sigma = 10)
#'
#' (x <- re(groups = 15, size = 15, replace = FALSE))
#' re(x, sigma = 1)
#'
#' (x <- re(groups = 2, size = 5, q = 2))
#' re(x, sigma = 10)
#'
#' (x <- 1:10)
#' re(x, sigma = 10)
#'
#' @export
re <- function (x, sigma, groups, size = NULL, q = 1, replace = TRUE) {
  if (!is.null(size)) {
    if (is.numeric(groups) & length(groups) == 1) {
      groups <- seq_len(groups)
      ngroups <- length(groups)
    }
    output <- rep(sample(groups, size = size, replace = replace), q) %>%
      factor(levels = groups)
  } else {
    if (!is.factor(x)) x <- factor(x)
    groups <- levels(x)
    beta <- rnorm(length(groups), 0, sigma)
    names(beta) <- groups
    output <- as.numeric(beta[x])
  }
  return(output)
}

exp_cor <- function (d, phi) {
  exp(-d/phi)
}

#' @title Simulate a Gaussian process
#'
#' @description
#' \code{gp} Simulate a spatial Gaussian process given a certain covariance function.
#'
#' @details
#' details.
#'
#' @param coords A list of coordinates
#' @param cor.model A character or function indicating the covariance function that
#' Should be used to compute the variance-covariance matrix
#' @param cor.params A list of the parameters required by the \code{cor.model} function.
#' @param sigma2 variance of the Gaussian process
#' @param size numeric value to simulate the covariate \code{x}
#'
#' @return A vector of the realization of the Gaussian Process
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' (x <- gp(list(s1 = NA), size = 10))
#' (s1 <- x[[1]])
#' # Simulate and plot the realization of a Gaussian process
#' (y <- gp(list(s1), cor.model = "exp_cor", cor.params = list(phi = 0.05)))
#'
#' (x <- gp(list(s1 = NA, s2 = NA), size = 10))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' # Simulate and plot the realization of a Gaussian process
#' (y <- gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.05)))
#' plot(s1, s2, cex = y)
#' # Plot with ggplot
#' # ggplot(data.frame(s1, s2, y), aes(s1, s2, col = y)) +
#' #  geom_point(size = 3)
#'
#' (x <- gp(list(s1 = "none", s2 = NA), size = 10))
#' (s1 <- x[[1]])
#' (s2 <- s1)
#' # Simulate and plot the realization of a Gaussian process
#' (y <- gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.05)))
#'
#' @importFrom stats dist rnorm runif
#' @importFrom purrr map
#'
#' @export
gp <- function (coords, cor.model, cor.params, sigma2 = 1, size = NULL) {
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(stats::runif(size)))
  } else {
    coords <- do.call(cbind, coords)
    n <- nrow(coords)
    distance <- as.matrix(dist(coords))
    varcov <- sigma2 * do.call(cor.model, c(list(distance), cor.params))
    right <- chol(varcov)
    output <- as.numeric(crossprod(right, rnorm(n)))
  }
  return(output)
}

