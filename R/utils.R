#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
#' @title Load system file for this package
#'
#' @description Load and return file from within the inst folder of this
#'   package.
#'
#' @param name name of file
#'
#' @export

malariaEqVivax_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'malariaEquilibriumVivax', mustWork = TRUE)
  ret <- readRDS(name_full)
}

