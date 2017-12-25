

#' @title Factor effects
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

#' @title Gaussian process effect
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
#' (x <- gp(s1 = NA, s2 = NA, beta = -1:1, size = 10))
#' (s1 <- x[[1]])
#' (s2 <- x[[2]])
#' gp(s1, s2, beta = -1:0)
#'
#' (s1 <- runif(10))
#' (x <- gp(s1, s2 = NA, beta = -1:1, size = 10))
#' (s2 <- x[[1]])
#' gp(s1, s2, beta = -1:0)
#'
#' @export
gp <- function (..., beta, size = NULL) {
  vars <- list(...)
  nvar <- purrr::map(vars, is.na) %>% do.call(sum, .)
  if (!is.null(size)) {
    output <- replicate(nvar, list(runif(size)))
  } else {
    output <- s1 * beta[1] + s2 * beta[2]
  }
  return(output)
}

