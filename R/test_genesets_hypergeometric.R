
#' Geneset ORA using hypergeometric test
#'
#' @description In most cases, it's more convenient to call the more generic `test_genesets` function which also deals with multiple-testing correction (per geneset source)
#'
#' It is assumed that the `genesets` and `genelist` parameters are in sync, i.e. `genesets` provided
#' here is the result of the `filter_genesets()` function (using the same `genelist` table)
#'
#' Only genesets with at least 1 significant gene are subjected to statistical testing (e.g. NA is returned for genesets without significant genes)
#'
#' @param genesets tibble with genesets, must contain columns 'id', 'ngenes' and 'ngenes_signif'
#' @param genelist tibble with genes, must contain column 'signif'. The number of rows in this table (where signif is not NA)
#' is assumed to be the total number of tested genes, the number of rows where signif==TRUE is assumed the total number of significant genes.
#' @param require_nsignif minimum number of 'signif genes' that overlap with a geneset; `NA` pvalues are returned for genesets with `ngenes_signif <= require_nsignif`.
#' This function 'prefilters' genesets, so beware that this will influence downstream multiple testing correction. Default is 1
#' @return input `genesets` table with results in the "pvalue" column
#' @seealso `test_genesets`
#' @export
test_genesets_hypergeometric = function(genesets, genelist, require_nsignif = 1L) {
  ngenes_signif = ngenes = pvalue = NULL # fix invisible bindings R package NOTE
  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns ngenes (integer type), ngenes_signif (integer type)" =
              length(genesets) > 0 && is.data.frame(genesets) &&  all(c("ngenes", "ngenes_signif") %in% colnames(genesets)) &&
              is.integer(genesets$ngenes) && is.integer(genesets$ngenes_signif))
  stopifnot("genelist parameter must be a non-empty data.frame/tibble with column 'signif' (logical/boolean type)" =
              length(genelist) > 0 && is.data.frame(genelist) && "signif" %in% colnames(genelist) &&
              is.logical(genelist$signif))
  stopifnot("require_nsignif parameter should be a single positive integer" = length(require_nsignif) == 1 && is.numeric(require_nsignif) && is.finite(require_nsignif) && require_nsignif >= 0)

  k_all = sum(!is.na(genelist$signif))     # total number of genes in user's genelist (background)
  k_signif = sum(genelist$signif %in% TRUE)  # total number of significant genes in user's genelist (foreground)

  genesets |>
    ungroup() |> # force ungrouping
    mutate(
      # for each geneset G, compute p-value from hypergeometric test
      pvalue = stats::phyper(
        ngenes_signif - 1,   # significant genes in geneset G - 1
        k_signif,            # all significant genes in test universe
        k_all - k_signif,    # genes in test universe that are not significant
        ngenes,            # all genes in geneset G
        lower.tail = FALSE),
      # set NA pvalue if geneset doesn't have sufficient overlap with significant genelist
      pvalue = ifelse(ngenes_signif < require_nsignif, NA, pvalue)
    )
}
