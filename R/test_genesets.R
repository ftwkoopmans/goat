
#' Perform geneset enrichment testing using any supported method
#'
#' @details
#' After application of the enrichment testing algorithm (e.g. GOAT, ORA or GSEA), multiple testing correction is applied to obtain adjusted p-values using `padjust_genesets`.
#' That function will first apply the specified pvalue adjustment procedure in the `padj_method` parameter within each 'source' in the genesets table. Second, it applies Bonferroni adjustment to all p-values according to the number of different geneset sources that were tested (or set `padj_sources = FALSE` to disable).
#'
#' For example, if the input is a genesets table that contains GO_CC, GO_BP and GO_MF genesets, first multiple testing correction is applied within each source (e.g. using FDR if so desired) and afterwards a Bonferroni correction is applied based on 3 repeated analyses.
#'
#' Note that this is more rigorous than typical GO tools; hypothetically, one could split all GO_CC pathways into 1000 different databases/'sources' and then run enrichment testing. Consequently, the multiple testing burden is reduced if one doesn't adjust p-values for the number of 'sources' as we do here.
#'
#' @param genesets tibble with genesets, must contain columns 'source', 'source_version', 'id', 'name', 'genes', 'ngenes', 'ngenes_signif'
#' @param genelist tibble with genes, must contain column 'gene' and 'test'. gene = character column, which are matched against list column 'genes' in genesets tibble. test = boolean column (you can set all to FALSE if not performing Fisher-exact or hypergeometric test downstream)
#' @param method method for overrepresentation analysis. Options: "goat", "hypergeometric", "fisherexact", "fisherexact_ease", "gsea"
#' @param padj_method first step of multiple testing correction; method for p-value adjustment, passed to `stats::p.adjust()` via `padjust_genesets()`, e.g. set "BH" to compute FDR adjusted p-values (default) or "bonferroni" for a more stringent procedure
#' @param padj_sources second step of multiple testing correction; apply Bonferroni adjustment to all p-values according to the number of geneset sources that were tested. Boolean parameter, set TRUE to enable (default) or FALSE to disable
#' @param padj_cutoff cutoff for adjusted p-value, `signif` column is set to TRUE for all values lesser-equals
#' @param padj_min_signifgenes if a value larger than zero is provided, this will perform additional post-hoc filtering; after p-value adjustment, set the `pvalue_adjust` to `NA` and `signif` to `FALSE` for all genesets with fewer than `padj_min_signifgenes` 'input genes that were significant' (`ngenes_signif` column in genesets table).
#' So this does not affect the accuracy of estimated p-values, in contrast to prefiltering genesets prior to p-value computation or adjusting p-values
#' @param ... further parameters are passed to the respective stats method
#' @return the input `genesets`, with results stored in columns 'pvalue', 'pvalue_adjust' and 'signif'
#' @export
test_genesets = function(genesets, genelist, method, padj_method = "BH", padj_sources = TRUE, padj_cutoff = 0.01, padj_min_signifgenes = 0L, ...) {
  algorithm = ngenes_signif = pvalue_adjust = pvalue = NULL # fix invisible bindings R package NOTE
  genesets = validate_genesets(genesets)
  genelist = validate_genelist(genelist)

  stopifnot("method parameter should be a single string, specifying one of the supported functions (see function documentation)" = length(method) == 1 && is.character(method))
  stopifnot("padj_method parameter should be a valid option as listed in `stats::p.adjust.methods`" = length(padj_method) == 1 && is.character(padj_method) && padj_method %in% stats::p.adjust.methods)
  stopifnot("padj_sources parameter should be a single boolean value" = length(padj_sources) == 1 && padj_sources %in% c(TRUE, FALSE))
  stopifnot("padj_cutoff parameter should be a single positive numeric value" = length(padj_cutoff) == 1 && is.numeric(padj_cutoff) && is.finite(padj_cutoff) && padj_cutoff >= 0)
  stopifnot("padj_min_signifgenes parameter should be a single positive integer (0 or higher)" = length(padj_min_signifgenes) == 1 && is.numeric(padj_min_signifgenes) && is.finite(padj_min_signifgenes) && padj_min_signifgenes >= 0)
  padj_min_signifgenes = as.integer(padj_min_signifgenes) # if supplied as 'numeric', round down and store as integer type
  result = NULL

  if(method == "goat" || method == "goat_precomputed") {
    result = test_genesets_goat_precomputed(genesets, genelist, ...)
  }
  if(method == "goat_fitfunction") {
    result = test_genesets_goat_fitfunction(genesets, genelist, ...)
  }
  if(method == "goat_bootstrap") {
    result = test_genesets_goat_bootstrap(genesets, genelist, ...)
  }
  if(method == "fisherexact") {
    result = test_genesets_fisherexact(genesets, genelist, ...)
  }
  if(method == "fisherexact_ease") {
    result = test_genesets_fisherexact(genesets, genelist, use_ease = T, ...)
  }
  if(method == "hypergeometric") {
    result = test_genesets_hypergeometric(genesets, genelist, ...)
  }
  if(method == "gsea") {
    result = test_genesets_gsea(genesets, genelist, ...)
  }

  # # fallthrough; guess if the user provided some custom function
  # if(is.null(result)) {
  #   f = tryCatch(match.fun(method, descend = FALSE), error = function(...) NULL)
  #   if(!is.function(f)) {
  #     if(grepl("::", method, fixed = T)) {
  #       f = tryCatch(utils::getFromNamespace(gsub(".*::", "", method), gsub("::.*", "", method), envir = .GlobalEnv), error = function(...) NULL)
  #     } else {
  #       f = tryCatch(utils::getFromNamespace(algorithm, "goat", envir = .GlobalEnv), error = function(...) NULL)
  #     }
  #   }
  #   if(is.function(f)) {
  #     result = f(genesets = genesets, genelist = genelist, ...)
  #     stopifnot("custom geneset-test-function should return the input 'genesets' table with an additional 'pvalue' numeric column included" =
  #                 is.data.frame(result) && "pvalue" %in% colnames(result) && is.numeric(result$pvalue))
  #   }
  # }


  # fallthrough check for valid method parameter
  stopifnot("unknown 'method' parameter" = !is.null(result))


  # finally, adjust p-values and flag significant
  result = result |>
    padjust_genesets(method = padj_method, cutoff = padj_cutoff, correct_sources = padj_sources) |>
    mutate(
      signif = signif & ngenes_signif >= padj_min_signifgenes,
      pvalue_adjust = ifelse(ngenes_signif < padj_min_signifgenes, NA, pvalue_adjust)
    ) |>
    score_geneset_oddsratio(genelist = genelist) |>
    arrange(pvalue)

  # settings as string
  settings... = parameters_prettyprint_length1(...)
  settings = sprintf("test_genesets(method='%s', padj_method='%s', padj_sources=%s, padj_cutoff=%s, padj_min_signifgenes=%s%s)",
                     method, padj_method, padj_sources, padj_cutoff, padj_min_signifgenes, ifelse(settings... == "", "", paste0(", ", settings...)))
  attr(result, "settings") <- c(attr(genesets, "settings"), settings)
  return(result)
}
