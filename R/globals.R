# R/globals.R
#' Declare data frame columns as global variables for R CMD check
#' @noRd
utils::globalVariables(c("ABC", "value", "Line", "x", "y", "xend", "yend"))
