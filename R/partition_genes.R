
#' Classify genes into 2 groups, e.g. to define significant or topN genes, resulting in a 'signif' column with boolean values
#'
#' This can be convenient to prepare the significant/test/foreground set for classical ORA,
#' e.g. `test_genesets()` with parameter `method = "fisherexact"`.
#'
#' @examples \dontrun{
#'   # significant hits
#'   genelist = partition_genes(genelist, col="pvalue_adjust", decreasing=F, cutoff=0.01)
#'   # abs(effectsize) >= 4
#'   genelist = partition_genes(genelist, col="effectsize", decreasing=T, use_abs=T, cutoff=4)
#'   # top 10% 'best' p-values
#'   genelist = partition_genes(genelist, col="pvalue", decreasing=F, fraction = 0.1)
#' }
#' @param genes gene tibble where each row is a unique gene, must contain column name `col`
#' @param col column name in `genes`
#' @param decreasing order `col` in descending (set TRUE) or ascending order (set FALSE, default) prior to partitioning?
#' @param use_abs use absolute values (default FALSE), e.g. when setting a threshold on effect-sizes
#' @param cutoff threshold for values in `col` to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @param fraction fraction of rows in `genes` tibble to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @param topn number of rows in `genes` tibble to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @export
partition_genes = function(genes, col, decreasing = FALSE, use_abs = FALSE, cutoff = NULL, fraction = NULL, topn = NULL) {
  signif = NULL # fix invisible bindings R package NOTE
  stopifnot("col parameter must be a string of length 1 (that matches a column in the genes table)" = length(col) == 1 && !is.na(col) && is.character(col) && col != "")
  stopifnot("decreasing parameter must be TRUE or FALSE" = length(decreasing) == 1 && decreasing %in% c(TRUE, FALSE))
  stopifnot("use_abs parameter must be TRUE or FALSE" = length(use_abs) == 1 && use_abs %in% c(TRUE, FALSE))
  stopifnot("cutoff parameter must be NULL or a positive numeric value" = is.null(cutoff) || (length(cutoff) == 1 && is.numeric(cutoff) && !is.na(cutoff) && cutoff > 0))
  stopifnot("fraction parameter must be NULL or a positive numeric value" = is.null(fraction) || (length(fraction) == 1 && is.numeric(fraction) && !is.na(fraction) && fraction > 0))
  stopifnot("topn parameter must be NULL or a positive numeric value" = is.null(topn) || (length(topn) == 1 && is.numeric(topn) && !is.na(topn) && topn > 0))

  stopifnot("genes parameter should be a nonempty data.frame that contains the column referred to in parameter col" =
              is.data.frame(genes) && nrow(genes) > 0 && col %in% colnames(genes) && is.numeric(genes |> pull(col)) )
  stopifnot("must specify 1, and only 1, of these parameters; cutoff, fraction, topn" = length(cutoff) + length(fraction) + length(topn) == 1)

  # - force ungrouping
  # - we need a temporary column to store the values to sort @ column `col`, which may have to be converted to absolute values
  # we'll use the `signif` column as that's our return value later, so this column name is available for overwriting
  genes = genes |> ungroup() |> mutate(signif = !!sym(col))
  if(use_abs) {
    genes = genes |> mutate(signif = abs(signif))
  }

  # value thresholding
  if(length(cutoff) == 1) {
    if(decreasing) { # decreasing = largest values on top = select all values larger-equals cutoff
      return(genes |> mutate(signif = signif >= cutoff))
    } else {
      return(genes |> mutate(signif = signif <= cutoff)) # sort ascending = take values lesser-equals cutoff
    }
  }

  # if we haven't returned at this point, user requested `fraction` or `topn`  -->>  both require sorting the table
  if(decreasing) {
    genes = genes |> arrange(desc(signif))
  } else {
    genes = genes |> arrange(signif)
  }

  if(length(fraction) == 1) {
    return(genes |> mutate(signif = seq_len(n()) <= n() * fraction))
  }

  if(length(topn) == 1) {
    return(genes |> mutate(signif = seq_len(n()) <= topn))
  }
}
