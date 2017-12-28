
#' @title Multivariate factor effects
#'
#' @description
#' \code{function} description.
#'
#' @details
#' \code{mfa} has two behaviours. If \code{size == NULL}, it evaluate the effects
#' applied to the argument \code{x}; otherwise, it generates a factor \code{x}.
#'
#' @param x factor vector to evaluate the effect
#' @param beta matrix of the effects for each level of \code{x}
#' @param labels character vector to name the \code{levels} of \code{x}
#' @param size numeric value to simulate the covariate \code{x}
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' (x <- mfa(beta = cbind(1:2, 2:3, 0:1), size = 10))
#' mfa(x, beta = cbind(1:2, 2:3, 0:1))
#'
#' (x <- mfa(beta = replicate(3, 0:2), size = 10))
#' mfa(x, beta = replicate(3, 0:2))
#'
#' (x <- mfa(beta = cbind(1, 2, 0), size = 10))
#' mfa(x, beta = cbind(1, 2, 0))
#'
#' @export
mfa <- function (x, beta, labels = 1:nrow(beta), size = NULL) {
  if (!is.null(size)) {
    output <- factor(sample(labels, size = size, replace = TRUE), levels = labels)
  } else {
    if (length(labels) > 1) {
      design <- model.matrix(~ -1 + x, data.frame(x))
    } else {
      design <- rep(1, length(x))
    }
  output <- as.numeric(design %*% beta)
  }
  return(output)
}

#' @title Multivariate Fixed Effect
#'
#' @description
#' \code{mfe} compute the multivariate fixed effect. Generally used with
#' \code{sim_model}.
#'
#' @details
#' details.
#'
#' @param x A vector of length \code{n} for which the fixed effect will be evaluated.
#' @param beta A vector of length \code{q}, this is the fixed effect for each response
#' variable.
#'
#' @return A matrix of dimension n x q.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' mfe(x = rnorm(10), beta = c(0.1, 0, 1))
#'
#' mfe(x = rnorm(10), beta = rep(10, 3))
#'
#' @export
mfe <- function (x, beta) {
  as.numeric(matrix(x) %*% matrix(beta, nrow = 1))
}

#' @title Multivariate random effects
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param x factor vector to evaluate the random effect
#' @param sigma variance matrix of the random effect
#' @param groups character vector to name the \code{levels} of \code{x} or a numeric
#' value indicating the number of groups of \code{x}
#' @param size numeric value to simulate the covariate \code{x}
#' @param replace An optional argument to simulate an independent random effect with
#' one repetition (e.g. \code{mre(groups = 100, size = 100, replace = FALSE)})
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mre(groups = 300, size = 300))
#' (effect <- mre(x, sigma = Sigma))
#'
#' id_uni <- match(unique(x), x)
#' effect_mat <- matrix(effect, ncol = nrow(Sigma))[id_uni, ]
#' cov(effect_mat)
#'
#' (x <- mre(groups = 300, size = 3000))
#' (effect <- mre(x, sigma = matrix(c(2, 1.5, 1.5, 2), nrow = 2)))
#'
#' id_uni <- match(unique(x), x)
#' effect_mat <- matrix(effect, ncol = 2)[id_uni, ]
#' cov(effect_mat)
#'
#' (x <- mre(groups = 10, size = 10, replace = FALSE))
#' (effect <- mre(x, sigma = matrix(c(2, 1.5, 1.5, 2), nrow = 2)))
#'
#' id_uni <- match(unique(x), x)
#' effect_mat <- matrix(effect, ncol = 2)[id_uni, ]
#' cov(effect_mat)
#'
#' @export
mre <- function (x, sigma, groups, size = NULL, replace = TRUE) {
  if (!is.null(size)) {
    if (is.numeric(groups) & length(groups) == 1) {
      groups <- seq_len(groups)
      ngroups <- length(groups)
    }
    output <- factor(sample(groups, size = size, replace = replace), levels = groups)
  } else {
    groups <- levels(x)
    right <- kronecker(chol(sigma), diag(length(groups)))
    beta <- crossprod(right, rnorm(length(groups) * nrow(sigma)))
    beta <- matrix(beta, nrow = length(groups))

    design <- model.matrix(~ -1 + x, data.frame(x))
    output <- as.numeric(design %*% beta)
  }
  return(output)
}


#' @title Simulate a Multivariate Gaussian process
#'
#' @description
#' \code{mgp} Simulate a Multivariate spatial Gaussian process known as linear model
#' of coregionalization Y(h) = AS(h), where S(h) is a vector of q independent Gaussian
#' processes.
#'
#' @details
#' details.
#'
#' @param coords A list of coordinates
#' @param variance A qxq cross-covariance matrix.
#' @param cor.model A character or function indicating the covariance function that
#' should be used to compute the correlation matrix
#' @param cor.params List of lists indicating the parameters for each response.
#' @param size numeric value to simulate the covariate \code{x}
#'
#' @return A vector of the realization of the Gaussian Process
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' # Generate coordinates
#' N <- 100
#' s1 <- 2 * runif(N)
#' s2 <- 2 * runif(N)
#'
#' # Covariance parameters
#' q <- 2
#' var <- sqrt(diag(c(4, 4)))
#' A <- matrix(c(1, - 0.8, 0, 0.6), nrow = 2)
#' variance <- var %*% tcrossprod(A) %*% var
#' cor.params <- list(list(phi = 0.08), list(phi = 0.08))
#'
#' # Generate the multivariate Gaussian process
#' y <- mgp(list(s1, s2), variance, "exp_cor", cor.params)
#' y1 <- y[1:N]
#' y2 <- y[(N + 1):(2 * N)]
#'
#' # Check correlation
#' cor(y1, y2)
#' plot(y1, y2)
#'
#' # Visualize the spatial
#' plot(s1, s2, cex = y1, col = 2)
#' points(s1, s2, cex = y2, col = 3)
#'
#' (x <- mgp(list(s1 = NA, s2 = NA), size = 100))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07))
#' variance = matrix(c(2, 1.5, 1.5, 2), nrow = 2)
#' (out <- mgp(list(s1, s2), variance, cor.model, cor.params))
#' out_mat <- matrix(out, ncol = 2)
#' plot(out_mat)
#' cov(out_mat)
#'
#' (x <- mgp(list(s1 = NA), size = 10))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' mgp(list(s1), variance = variance, cor.model = cor.model, cor.params = cor.params)
#'
#' (x <- mgp(list(s1 = NA), size = 100))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' mgp(list(s1), variance, cor.model, cor.params)
#'
#' @importFrom stats dist runif rnorm model.matrix
#' @importFrom purrr map reduce
#'
#' @export

mgp <- function (coords, variance, cor.model, cor.params, size = NULL) {
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(stats::runif(size)))
  } else {
    coords <- do.call(cbind, coords)
    n <- nrow(coords)
    q <- nrow(variance)
    distance <- as.matrix(dist(coords))
    rights <- purrr::map(cor.params, ~ chol(do.call(cor.model, c(list(distance), .))))

    igp <- purrr::map(rights, ~ as.numeric(crossprod(., stats::rnorm(n)))) %>%
      purrr::reduce(rbind) %>%
      as.data.frame() %>%
      as.list()

    variance_chol <- chol(variance)
    mgp <- purrr::map(igp, ~ as.numeric(crossprod(variance_chol, .))) %>%
      purrr::reduce(rbind)
    output <- as.numeric(mgp)
  }
  return(output)
}
