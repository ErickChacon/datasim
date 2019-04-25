
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
#' \code{mgp} handles the evaluation of multivariate Gaussian processes with two
#' behaviours. It evaluates the multivariate effects applied to the coordinates
#' \code{coords} if \code{size == NULL} or it simulates the coordinates
#' \code{coords} if \code{size} is provided.
#'
#' @details
#' The \eqn{q}-dimensional multivariate Gaussian process \eqn{Y(h)} is
#' expressed as a linear combination
#' \deqn{Y(h) = AS(h).}
#' where \eqn{S(h)} is a vector of \eqn{m} standardized independent Gaussian
#' processes at location \eqn{h} and \eqn{A} a \eqn{qxm} matrix of coefficients.
#'
#' @param coords A list of replicated coordinates. Each coordinate is a \eqn{q} replicated
#' vector of the unique coordinates. If \code{size != NULL}, \code{coords} is the
#' output of the function \code{mgp}.
#' @param A A qxm matrix that defines the relationship between \eqn{Y(h)} and
#' \eqn{S(h)}.
#' @param cor.model A character or function indicating the correlation function to
#' be used to compute the correlation matrix
#' @param cor.params A nested list of the parameters of the correlation function for
#' each response.
#' @param size A numeric value \eqn{n} representing the number of locations, it is used to
#' simulate the coordinates \code{coords}. In case \code{size == NULL}, \code{mgp}
#' simulates the multivariate Gaussian process.
#' @param range Range of the coordinates of the Gaussian processes.
#'
#' @return A list of simulated replicated coordinates in case \code{size} is provided;
#' otherwise, a \eqn{nq}-length numeric vector of the evaluated multivariate Gaussian
#' process.
#'
#' @author Erick A. Chacón-Montalván
#'
#' @examples
#' ## Simulation of a temporal multivariate Gaussian process
#'
#' # Define the number of location and responses, and a linear transformation matrix
#' n <- 200
#' q <- 3
#' A <- matrix(c(1, 0, 1, 0, 1, - 0.8), nrow = q)
#'
#' # Simulate coordinates
#' coords <- mgp(list(time = NA), A, size = n)
#'
#' # Simulate a multivariate Gaussian process
#' cor.params <- list(list(phi = 0.08), list(phi = 0.15))
#' y <- mgp(coords, A, "exp_cor", cor.params)
#'
#' # Visualize the temporal multivariate Gaussian process
#' data <- data.frame(time = coords[[1]], y, response = factor(rep(1:q, each = n)))
#' ggplot(data, aes(time)) +
#'   geom_line(aes(y = y, col = response))
#'
#' ## Simulation of a spatial multivariate Gaussian process
#'
#' # Define the number of location and responses, and a linear transformation matrix
#' n <- 1500
#' q <- 2
#' A <- matrix(c(1, -1), nrow = q)
#'
#' # Simulate coordinates
#' coords <- mgp(list(s1 = NA, s2 = NA), A, size = n)
#'
#' # Simulate a multivariate Gaussian process
#' cor.params <- list(list(phi = 0.1))
#' y <- mgp(coords, A, "exp_cor", cor.params)
#'
#' # Visualize the temporal multivariate Gaussian process
#' data <- data.frame(s1 = coords[[1]], s2 = coords[[2]], y,
#'   response = factor(rep(1:q, each = n)))
#' ggplot(data, aes(s1, s2)) +
#'   geom_point(aes(size = y, col = y)) +
#'   facet_wrap(~ response)
#'
#' @importFrom stats dist runif rnorm model.matrix
#' @importFrom purrr map reduce
#'
#' @export
mgp <- function (coords, A, cor.model, cor.params, size = NULL, range = 1) {
  q <- nrow(A)
  if (!is.null(size)) {
    ncoords <- purrr::map(coords, is.na) %>% do.call(sum, .)
    output <- replicate(ncoords, list(rep(range * stats::runif(size), q)))
  } else {
    coords <- do.call(cbind, coords)
    coords <- coords[1:(nrow(coords) / q), , drop = FALSE]
    n <- nrow(coords)
    distance <- as.matrix(dist(coords))
    rights <- purrr::map(cor.params, ~ chol(do.call(cor.model, c(list(distance), .))))

    igp <- purrr::map(rights, ~ as.numeric(crossprod(., stats::rnorm(n)))) %>%
      do.call(cbind, .)
    mgp <- tcrossprod(igp, A)
    output <- as.numeric(mgp)
  }
  return(output)
}

