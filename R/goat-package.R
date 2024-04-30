#' goat: Gene Set Analysis Using the Gene Set Ordinal Association Test
#'
#' @description
#' Perform gene set enrichment analyses using the
#' Gene set Ordinal Association Test (GOAT) algorithm
#' and visualize your results.
#' @name goat
#' @useDynLib goat, .registration=TRUE
#' @keywords internal
#' @import dplyr
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

#' cleanup Rcpp code
#'
#' @param libpath library path
#' @noRd
.onUnload = function(libpath = NULL) {
  library.dynam.unload("goat", libpath)
}
