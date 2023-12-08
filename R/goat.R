#' goat: Perform geneset enrichment and overrepresentation analyses using various algorithms, including GOAT.
#'
#' @description
#' TODO: package description
#'
#' @docType package
#' @name goat
#' @import dplyr
#' @importFrom Rcpp sourceCpp
#' @useDynLib goat, .registration=TRUE
NULL



#' cleanup Rcpp code
#'
#' @param libpath library path
.onUnload = function(libpath = NULL) {
  library.dynam.unload("goat", libpath)
}
