#' \code{facetsSuite} package
#'
#' Functions to run and parse output from \href{https://github.com/mskcc/facets}{FACETS}.
#'
#' @docType package
#' @name facetsSuite
#' @importFrom dplyr %>%
#' @import data.table
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."), add = TRUE)