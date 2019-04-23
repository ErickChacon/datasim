
#' @title Multivariate Factor Effects
#'
#' @description
#' \code{mfa} handles the evaluation of multivariate factor effects with two
#' behaviours. It evaluates the multivariate effects applied to the replicated factor
#' \code{x} if \code{size == NULL} or simulates the replicated factor \code{x} if
#' \code{size} is provided.
#'
#' @details
#' Considering \eqn{x*} the \eqn{n}-length factor under study with \eqn{k} levels and
#' \eqn{X} the associated \eqn{nxk} design matrix with dummy variables corresponding
#' to each level of the factor, the returning multivariate effect is
#' \deqn{vec(XB),}
#' where \eqn{B} is a \eqn{kxq} matrix of regression coefficients and \eqn{vec(.)}
#' represents a vectorization by columns of the provided matrix.
#'
#' @param x A replicated covariate factor of length \eqn{nq} to evaluate the
#' multivariate effects. It is a \eqn{q} times replicated vector of \eqn{x*}. If
#' \code{size != NULL}, \code{x} is the output of the function.
#' @param beta A \eqn{kxq} matrix B of regression coefficients representing the effects
#' for each level of \code{x} on the response variables. Each column represents the
#' effect to each response variable.
#' @param levels A character vector of length \eqn{k} to name the \code{levels} of the
#' factor \code{x}.
#' @param size A numeric value \eqn{n} representing the number of units, it is used to
#' simulate the covariate \code{x}. In case \code{size == NULL}, \code{mfa} evaluates
#' the effects.
#'
#' @return A simulated replicated factor \eqn{x} in case \code{size} is provided;
#' otherwise, a \eqn{nq}-length numeric vector of the evaluated multivariate
#' effects.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' # Different effects for each response.
#' (x <- mfa(beta = cbind(1:2, 2:3, 0:1), levels = c("F", "M"), size = 10))
#' mfa(x, beta = cbind(1:2, 2:3, 0:1))
#'
#' # Same effects for each response.
#' (x <- mfa(beta = replicate(3, 0:2), size = 10))
#' mfa(x, beta = replicate(3, 0:2))
#'
#' # Differrent intercepts for each response.
#' (x <- mfa(beta = cbind(1, 2, 0), size = 10))
#' mfa(x, beta = cbind(1, 2, 0))
#'
#' @export
mfa <- function (x, beta, levels = 1:nrow(beta), size = NULL) {
  q <- ncol(beta)
  if (!is.null(size)) {
    output <- rep(sample(levels, size = size, replace = TRUE), q) %>%
      factor(levels = levels)
  } else {
    if (!is.factor(x)) x <- factor(x)
    x <- x[1:(length(x) / q)]
    if (length(levels) > 1) {
      design <- model.matrix(~ -1 + x, data.frame(x))
    } else {
      design <- rep(1, length(x))
    }
  output <- as.numeric(design %*% beta)
  }
  return(output)
}

#' @title Multivariate Fixed Effects
#'
#' @description
#' \code{mfe} evaluates the multivariate fixed effects for a continuous replicated
#' covariate \code{x}.
#'
#' @details
#' Considering \eqn{x*} the continuous covariate of length \eqn{n}, the returning
#' multivariate effect is
#' \deqn{vec(x*b),}
#' where \eqn{b} is a \eqn{q}-length row vector of regression coefficients and
#' \eqn{vec(.)} represents a vectorization by columns of the provided matrix.
#'
#' @param x A replicated continuous covariate of length \code{nq} for which the
#' fixed effects are evaluated. It is a \eqn{q}-times replicated vector of \eqn{x*}.
#' @param beta A vector of length \code{q}, where each element is a regression
#' coefficient associated to each response variable.
#'
#' @return A \eqn{nq}-length numeric vector of the evaluated multivariate effects.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' # Different effects for each response.
#' mfe(x = rep(rnorm(5), 3), beta = c(0.1, 0, 1))
#'
#' # Same effects for each response.
#' mfe(x = rep(rnorm(5), 3), beta = rep(10, 3))
#'
#' @export
mfe <- function (x, beta) {
  np <- length(x)
  p <- length(beta)
  n <- np / p
  output <- x * rep(beta, each = n)
  # as.numeric(matrix(x) %*% matrix(beta, nrow = 1))
  return(output)
}

#' @title Multivariate Random Effects
#'
#' @description
#' \code{mre} handles the evaluation of multivariate random effects with two
#' behaviours. It evaluates the multivariate effects applied to the factor \code{x}
#' if \code{size == NULL} and it simulates the replicated factor \code{x} if
#' \code{size} is provided.
#'
#' @details
#' Considering \eqn{x*} the \eqn{n}-length factor under study with \eqn{k} levels and
#' \eqn{X} the associated \eqn{nxk} design matrix with dummy variables corresponding
#' to each level of the factor, the returning multivariate effect is
#' \deqn{vec(XU),}
#' where \eqn{U} is a \eqn{kxq} matrix of random effects and \eqn{vec(.)}
#' represents a vectorization by columns of the provided matrix. Each row \eqn{u_i}
#' of \eqn{U} is assumed to come from a zero-mean normal distribution with covariance
#' matrix \eqn{S} of dimension \eqn{qxq},
#'
#' @param x A replicated covariate factor of length \eqn{nq} to evaluate the
#' multivariate effects. It is a \eqn{q} times replicated vector of \eqn{x*}. If
#' \code{size != NULL}, \code{x} is the output of the function.
#' @param sigma A \eqn{qxq} covariance matrix \eqn{S} for the random effects.
#' @param groups A character vector of length \eqn{k} to name the \code{levels} of
#' the factor \code{x} or a numeric value indicating the number of groups of
#' \code{x}.
#' @param size A numeric value \eqn{n} representing the number of units, it is used to
#' simulate the covariate \code{x}. In case \code{size == NULL}, \code{mfa} evaluates
#' the effects.
#' @param replace An logical value provided to the function \code{sample} to allow
#' repetition of groups or not. It is used to simulate an independent random effect with
#' one repetition (e.g. \code{mre(groups = 100, size = 100, replace = FALSE)}).
#'
#' @return A simulated replicated factor \eqn{x} in case \code{size} is provided;
#' otherwise, a \eqn{nq}-length numeric vector of the evaluated multivariate
#' effects.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' # Multivariate random effects for 10 units belonging to 3 groups.
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mre(groups = 3, sigma = Sigma, size = 10))
#' (effect <- mre(x, sigma = Sigma))
#'
#' # Multivariate independent random effects for 10 units.
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mre(groups = 10, sigma = Sigma, size = 10, replace = FALSE))
#' (effect <- mre(x, sigma = Sigma))
#'
#' # Multivariate random effects for 500 units belonging to 100 groups.
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mre(groups = 100, sigma = Sigma, size = 500))
#' (effect <- mre(x, sigma = Sigma))
#'
#' # Check empirical covariance matrix of the random effects.
#' cov(unique(matrix(effect, ncol = nrow(Sigma))))
#'
#' @export
mre <- function (x, sigma, groups, size = NULL, replace = TRUE) {
  q <- nrow(sigma)
  if (!is.null(size)) {
    if (is.numeric(groups) & length(groups) == 1) {
      groups <- seq_len(groups)
      ngroups <- length(groups)
    }
    output <- rep(sample(groups, size = size, replace = replace), q) %>%
      factor(levels = groups)
  } else {
    x <- x[1:(length(x) / q)]
    if (!is.factor(x)) x <- factor(x)
    groups <- levels(x)
    right <- kronecker(chol(sigma), diag(length(groups)))
    beta <- crossprod(right, rnorm(length(groups) * nrow(sigma)))
    beta <- matrix(beta, nrow = length(groups))

    design <- model.matrix(~ -1 + x, data.frame(x))
    output <- as.numeric(design %*% beta)
  }
  return(output)
}

#' @title Multivariate Gaussian Processes
#'
#' @description
#' \code{mgp} Simulate a Multivariate spatial Gaussian process known as linear model
#' of coregionalization Y(h) = AS(h), where S(h) is a vector of q independent Gaussian
#' processes.
#'
#' @details
#' details.
#'
#' @param coords A list of coordinates.
#' @param variance A qxq cross-covariance matrix Sigma such as \code{A=chol(Sigma)}.
#' @param cor.model A character or function indicating the covariance function that
#' should be used to compute the correlation matrix
#' @param cor.params List of lists indicating the parameters for each response.
#' @param size numeric value to simulate the covariate \code{x}
#' @param range range on the coordinates of the Gaussian process
#'
#' @return A vector of the realization of the Gaussian Process
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#'
#' # Generate coordinates
#' N <- 100
#' s1 <- rep(2 * runif(N), 2)
#' s2 <- rep(2 * runif(N), 2)
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
#' variance = matrix(c(2, 1.5, 1.5, 2), nrow = 2)
#' (x <- mgp(list(s1 = NA, s2 = NA), variance, size = 100))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07))
#' (out <- mgp(list(s1, s2), variance, cor.model, cor.params))
#' out_mat <- matrix(out, ncol = 2)
#' plot(out_mat)
#' cov(out_mat)
#'
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mgp(list(s1 = NA), variance, size = 10))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' mgp(list(s1), variance = variance, cor.model = cor.model, cor.params = cor.params)
#'
#' variance <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' (x <- mgp(list(s1 = NA), variance, size = 100))
#' (s1 <- x[[1]])
#' cor.model <- "exp_cor"
#' cor.params <- list(list(phi = 0.05), list(phi = 0.07), list(phi = 1))
#' mgp(list(s1), variance, cor.model, cor.params)
#'
#' @importFrom stats dist runif rnorm model.matrix
#' @importFrom purrr map reduce
#'
#' @export

mgp <- function (coords, variance, cor.model, cor.params, size = NULL, range = 1) {
  q <- nrow(variance)
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(rep(range * stats::runif(size), q)))
  } else {
    coords <- do.call(cbind, coords)
    coords <- coords[1:(nrow(coords) / q), , drop = FALSE]
    n <- nrow(coords)
    # q <- nrow(variance)
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
