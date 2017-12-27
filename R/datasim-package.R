
#' @title datasim: Data simulation based on models formulas
#'
#' @docType package
#' @name day2day
#'
#' @description
#' The \code{datasim} package provide tools to simulate data frames based on models
#' formulas. It works with different types of effects and multivariate models.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @importFrom magrittr %>%
"_PACKAGE"

if(getRversion() >= "2.15.1") utils::globalVariables(
  c(".",
    "type", "type_order", "covs"
    )
  )

