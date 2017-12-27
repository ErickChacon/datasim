
#' @title Multivariate factor effects
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
#' \code{msim_model}.
#'
#' @details
#' details.
#'
#' @param x A vector of length n for which the fixed effect will be evaluated.
#' @param beta A vector of length q, this is the fixed effect for each response
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
#' @param par.
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
#' @param s1 First coordinate
#'
#' @param s2 Second coordinate
#'
#' @param cov.model A character or function indicating the covariance function that
#' Should be used to compute the variance-covariance matrix
#'
#' @param variance A qxq matrix of non-spatial covariance.
#' @param nugget A qxq diagonal matrix of non-spatial noise.
#' @param phi A q-length vector of decay parameters.
#' @param kappa A q-length vector of kappa parameters if Matern spatial correlation
#' function is used.
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
#' nugget <- diag(0, q)
#' phi <- rep(1 / 0.08, q)
#'
#' # Generate the multivariate Gaussian process
#' y <- mgp(s1, s2, "exponential", variance, nugget, phi)
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
#'
#' (x <- mgp(s1 = NA, s2 = NA, size = 100))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' coords <- cbind(s1, s2)
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07))
#' variance = matrix(c(2, 1.5, 1.5, 2), nrow = 2)
#' (out <- mgp(s1, s2, variance = variance, cor.model = cor.model, cor.params = cor.params))
#' out_mat <- matrix(out, ncol = 2)
#' plot(out_mat)
#' cov(out_mat)
#'
#' (x <- mgp(s1 = NA, size = 10))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' mgp(s1, variance = variance, cor.model = cor.model, cor.params = cor.params)
#'
#' (x <- mgp(s1 = NA, size = 100))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' mgp(s1, variance = variance, cor.model = cor.model, cor.params = cor.params)
#'
#' @importFrom spBayes mkSpCov
#'
#' @export

mgp <- function (..., variance = NULL, cor.model = NULL, cor.params = NULL,
                 size = NULL) {
  coords <- list(...)
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(runif(size)))
  } else {
    coords <- do.call(cbind, coords)
    n <- nrow(coords)
    q <- nrow(variance)
    distance <- as.matrix(dist(coords))
    rights <- map(cor.params, ~ chol(do.call(cor.model, c(list(distance), .))))

    igp <- map(rights, ~ as.numeric(crossprod(., rnorm(n)))) %>%
      reduce(rbind) %>%
      as.data.frame() %>%
      as.list()

    variance_chol <- chol(variance)
    mgp <- map(igp, ~ as.numeric(crossprod(variance_chol, .))) %>% reduce(rbind)
    output <- as.numeric(mgp)
  }
  return(output)
}
