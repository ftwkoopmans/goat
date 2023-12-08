
#' Compute odds-ratio for each geneset
#'
#' gs_signif = number of significant genes in geneset G that intersect with user's genelist (i.e. foreground genes in G)
#' gs_all  = number of genes in geneset G that intersect with user's genelist (i.e. foreground+background genes in G)
#' k_signif = total number of _significant_ genes in user's genelist
#' k_all  = total number of genes in user's genelist
#'
#' gs_signif/gs_all = ratio of foreground genes in geneset G
#' k_signif/k_all = ratio of overall foreground genes (i.e. expected value for a random geneset)
#'
#' oddsratio = (gs_signif/gs_all) / (k_signif/k_all)
#'
#' @param genesets tibble with genesets, must contain columns 'source', 'id', 'genes'
#' @param genelist tibble with genes, must contain column 'gene', 'log2fc', 'pvalue'
#' @export
score_geneset_oddsratio = function(genesets, genelist) {
  ngenes_signif = ngenes = NULL # fix invisible bindings R package NOTE
  stopifnot(length(genesets) > 0 && is.data.frame(genesets) && all(c("ngenes", "ngenes_signif") %in% colnames(genesets)))
  stopifnot(length(genelist) > 0 && is.data.frame(genelist) && "signif" %in% colnames(genelist))

  k_all = nrow(genelist)                 # total number of genes in user's genelist
  k_signif = sum(genelist$signif %in% TRUE)  # total number of significant genes in user's genelist
  genesets |> mutate(score_oddsratio = (ngenes_signif/ngenes) / (k_signif/k_all))
}
