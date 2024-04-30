
#' Classify genes into 2 groups, e.g. to define significant or topN genes, resulting in a 'signif' column with boolean values
#'
#' This can be convenient to prepare the significant/test/foreground set for classical ORA,
#' e.g. `test_genesets()` with parameter `method = "fisherexact"`. Note that the GOAT geneset
#' enrichment algorithm does not use data in the 'signif' column of the input genelist.
#'
#' @examples \donttest{
#' # note: this example will download 1 files of approx 4MB
#'
#' # store the downloaded files in the following directory. Here, the temporary file
#' # directory is used. Alternatively, consider storing this data in a more permanent location.
#' # e.g. output_dir="~/data/goat" on unix systems or output_dir="C:/data/goat" on Windows
#' output_dir = tempdir()
#'
#' # Download an example gene list, i.e. one of the datasets analyzed in the GOAT manuscript.
#' datasets = download_goat_manuscript_data(output_dir)
#' genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`
#'
#' # example 1: significant hits
#' genelist = partition_genes(genelist, col="pvalue_adjust", decreasing=FALSE, cutoff=0.01)
#' cat(sum(genelist$signif), "/", nrow(genelist), "are signif\n")
#'
#' # example 2: abs(effectsize) >= 5
#' genelist = partition_genes(genelist, col="effectsize", decreasing=TRUE, use_abs=TRUE, cutoff=5)
#' cat(sum(genelist$signif), "/", nrow(genelist), "are signif\n")
#'
#' # example 3: top 10% 'best' p-values
#' genelist = partition_genes(genelist, col="pvalue", decreasing=FALSE, fraction = 0.1)
#' cat(sum(genelist$signif), "/", nrow(genelist), "are signif\n")
#' }
#' @param genes gene tibble where each row is a unique gene, must contain column name `col`
#' @param col column name in `genes`
#' @param decreasing order `col` in descending (set TRUE) or ascending order (set FALSE, default) prior to partitioning?
#' @param use_abs use absolute values (default FALSE), e.g. when setting a threshold on effect-sizes
#' @param cutoff threshold for values in `col` to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @param fraction fraction of rows in `genes` tibble to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @param topn number of rows in `genes` tibble to select  (must provide exactly 1 parameter for filtering, either `cutoff`, `fraction` or `topn`)
#' @return input table `genes` with results in the "signif" column
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
