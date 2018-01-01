
#' @title datasim: Data simulation based on models formulas
#'
#' @docType package
#' @name datasim
#'
#' @description
#' The \code{datasim} package provide tools to simulate data frames based on models
#' formulas. It works with different types of effects and multivariate models. The
#' user define a list of formulas than later is used to i) simulate the covariates
#' presented in the formula and to ii) simulate the response variable.
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


#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
