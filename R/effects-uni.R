
#' @title Univariate factor effects
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
#' @param par.
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
#' @export
re <- function (x, sigma, groups, size = NULL) {
  if (!is.null(size)) {
    if (is.numeric(groups) & length(groups) == 1) {
      groups <- seq_len(groups)
      ngroups <- length(groups)
    }
    output <- factor(sample(groups, size = size, replace = TRUE), levels = groups)
  } else {
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
#' @param s1 First coordinate
#'
#' @param s2 Second coordinate
#'
#' @param cov.model A character or function indicating the covariance function that
#' Should be used to compute the variance-covariance matrix
#'
#' @param cov.params A list of the parameters required by the \code{cov.model} function.
#'
#' @return A vector of the realization of the Gaussian Process
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' (x <- gp(s1 = NA, size = 10))
#' (s1 <- x[[1]])
#' # Simulate and plot the realization of a Gaussian process
#' (y <- gp(s1, cor.model = "exp_cor", cor.params = list(phi = 0.05)))
#'
#' (x <- gp(s1 = NA, s2 = NA, size = 10))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' # Simulate and plot the realization of a Gaussian process
#' (y <- gp(s1, s2, cor.model = "exp_cor", cor.params = list(phi = 0.05)))
#' plot(s1, s2, cex = y)
#' # Plot with ggplot
#' # ggplot(data.frame(s1, s2, y), aes(s1, s2, col = y)) +
#' #  geom_point(size = 3)
#'
#' @importFrom stats dist rnorm
#'
#' @export
gp <- function (..., cor.model = NULL, cor.params = NULL, sigma2 = 1, size = NULL) {
  coords <- list(...)
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(runif(size)))
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

