
#' filter a geneset table; intersect with an array of genes-of-interest then apply cutoffs on min/max genes per geneset
#'
#' @param genesets tibble with genesets, must contain columns 'id', 'genes' and 'ngenes'
#' @param genelist tibble with genes, must contain column 'gene' and 'signif'. gene = character column, which are matched against list column 'genes' in genesets tibble. signif = boolean column (you can set all to FALSE if not performing Fisher-exact or hypergeometric test downstream)
#' @param min_overlap integer, minimum number of genes in the `genelist` table that must match a geneset. Must be at least 1
#' @param max_overlap integer, maximum number of genes in the `genelist` table that must match a geneset. Set to NA to disable
#' @param max_overlap_fraction analogous to `max_overlap`, which limits the max geneset size to a given N, this parameter defines the maximum geneset size that is to be retained as a fraction of the input genelist length. For example, setting this to 0.5 will remove all genesets that contain more than half the genes in the input genelist (i.e. testing enrichment of a geneset that contains 1000 out of a total 1200 genes from your input genelist is probably meaningless). Defaults to 50%
#' @param min_signif integer, minimum number of genes in the `genelist` table that are `signif==TRUE` and match a geneset. Be careful, this is "prefiltering" and will affect the correcteness / calibration of estimated geneset p-values. For GOAT and GSEA, this is NOT RECOMMENDED. Set to NA to disable (default)
#' @param max_size integer, maximum number of genes in the geneset (i.e. prior to intersect with user's gene list provided as `genelist`). Optionally, use this to remove highly generic terms. Set to NA to disable
#' @param dedupe boolean, remove duplicate genesets (as determined after intersection with `genelist`)
#' @return the input `genesets` filtered for the subset of rows that match user's filter parameters
#' @export
filter_genesets = function(genesets, genelist, min_overlap = 10L, max_overlap = 1500L, max_overlap_fraction = 0.5, min_signif = NA, max_size = NA, dedupe = FALSE) {
  ngenes = genes = ngenes_input = isdupe = signif = gene = genes_signif = ngenes_signif = NULL # fix invisible bindings R package NOTE
  genesets = validate_genesets(genesets, require_signif = FALSE)
  genelist = validate_genelist(genelist)

  stopifnot(length(min_overlap) == 1 && ((is.numeric(min_overlap) & min_overlap > 0) || is.na(min_overlap)))
  stopifnot(length(max_overlap) == 1 && ((is.numeric(max_overlap) & max_overlap > 0) || is.na(max_overlap)))
  stopifnot(length(min_signif) == 1 && ((is.numeric(min_signif) & min_signif >= 0) || is.na(min_signif)))
  stopifnot(length(max_size) == 1 && ((is.numeric(max_size) & max_size > 0) || is.na(max_size)))
  if(!is.finite(min_overlap)) min_overlap = 1L   # users can provide NA to disable filtering
  if(!is.finite(max_overlap)) max_overlap = Inf # users can provide NA to disable filtering
  if(!is.finite(min_signif)) min_signif = 0L # users can provide NA to disable filtering
  if(!is.finite(max_size)) max_size = Inf # users can provide NA to disable filtering

  if(min_signif > 0) {
    cat("Warning: the 'min_signif' parameter is enabled. Be careful, this is \"prefiltering\" and will affect the correcteness / calibration of estimated geneset p-values. For GOAT and GSEA, this is NOT RECOMMENDED\n")
  }

  # settings as string
  settings = sprintf("filter_genesets(min_overlap=%s, max_overlap=%s, max_overlap_fraction=%s, min_signif=%s, max_size=%s, dedupe=%s)",
                     min_overlap, max_overlap, max_overlap_fraction, min_signif, max_size, dedupe)

  # optionally we can remove genesets that contain more than half of the genes in the input genelist
  # (i.e. does enrichment testing a geneset with 1000 genes in a genelist of 1200 genes make sense ?!)
  max_overlap = as.integer(min(max_overlap, ceiling(nrow(genelist) * max_overlap_fraction)))

  # l = list of arrays. e.g. list(letters[1:3], letters[2:4], letters[3:1])
  finddupes = function(l) {
    if(length(l) == 0) return()
    if(length(l) == 1) return(F)
    duplicated.default(lapply(l, sort))
  }

  # - force ungrouping
  # - up-front filtering of genesets that we won't use regardless of intersection with foreground
  x = genesets |> ungroup() |> filter(ngenes >= min_overlap)
  # optional filter
  if(is.finite(max_size)) {
    x = x |> filter(ngenes <= max_size)
  }

  x = x |>
    # fast intersection of background and foreground gene lists by unlisting, vectorized match, relisting
    tidyr::unnest(genes) |> # to long format
    filter(genes %in% genelist$gene) |> # retain only geneset*gene mappings matching the foreground set
    tidyr::chop(genes) |> # collapse / back to genes as a nested list
    # annotated gene count  __prior to intersection with user gene list__
    select(-tidyselect::any_of(c("ngenes_input", "ngenes_signif", "genes_signif"))) |>
    tibble::add_column(ngenes_input = 0, .before = "ngenes") |>
    # update gene count column, it now represents number of genes in geneset that intersect with user's gene list
    mutate(
      ngenes_input = ngenes,
      ngenes = lengths(genes)
    ) |>
    # filter by size (number of overlapping genes in genelist * geneset)
    filter(ngenes >= min_overlap & ngenes <= max_overlap) |>
    tibble::add_column(ngenes_signif = 0, .after = "ngenes") |>
    # this'll be a list type column, but don't have to initialize as one (e.g. vector("list", length(ngenes)) )
    tibble::add_column(genes_signif = 0, .after = "genes")

  # deduplication within each group of genesets with the same length
  if(dedupe) {
    x = x |>
      arrange(ngenes_input) |> # ascending, i.e. smallest geneset (originally) on top
      group_by(ngenes) |> # can only be a dupe if vector length is equal, so efficiently check within same-length
      mutate(isdupe = finddupes(genes)) |>
      ungroup() |>
      filter(isdupe == F) |>
      select(-isdupe)
  }

  # gene IDs that are to be tested in downstream hypergeometric or Fisher-exact tests
  gene_signif = genelist |> filter(signif %in% TRUE) |> pull(gene) # use %in% to deal with potential NA values
  x = x |>
    mutate(
      genes_signif = lapply(genes, intersect, gene_signif),
      ngenes_signif = lengths(genes_signif)
    )

  # optionally, an additional filter applied to the subset of genes in the universe that are to be used as foreground
  if(min_signif > 0) {
    x = x |> filter(ngenes_signif >= min_signif)
  }

  if(nrow(x) == 0) {
    cat("filter_genesets() yields an empty result !\nAre the gene identifiers in your 'genesets' and 'genelist tables of the same type? e.g. both tables should contain NCBI Entrez gene IDs, or both use HGNC identifiers, or Ensembl gene IDs. Another common mistake is using different species, so double-check that both tables contain e.g. human gene identifiers\n")
  }

  attr(x, "settings") <- c(attr(genesets, "settings"), settings)
  return(x)
}
