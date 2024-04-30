
#' GSEA as implemented in the fgsea R package
#'
#' @description In most cases, it's more convenient to call the more generic `test_genesets` function which also deals with multiple-testing correction (per geneset source)
#'
#' @param genesets data.frame/tibble with geneset and gene columns
#' @param genelist data.frame/tibble with gene and score columns. Should contain columns gene and either pvalue or effectsize, depending on `score_type` parameter
#' @param score_type how to compute gene scores? options: "pvalue", "effectsize", "custom".
#' Option "pvalue" uses -log10 transformed values from the pvalue column in `genelist`.
#' Option "effectsize" uses values from the effectsize column in `genelist` as-is.
#' Option "custom" expects 2 additional parameters; `gsea_genelist_col` should be a column name in `genelist` to be used for fGSEA (values used as-is), `gsea_scoretype` should be the respective value for the fGSEA scoreType parameter ('pos', 'neg' or 'std)
#' @param parallel_threads number of threads to use for parallel processing. Set to 0 to automatically select all available processors/cores, set to 1 to disable (default) or to N to use N processes.
#' Note that multiprocessing sometimes breaks on RStudio + Windows, hence this parameter is set to 1 to disable multiprocessing by default for now
#' @param gseaParam passed to `fgsea::fgsea()`, from manual: "GSEA parameter value, all gene-level statis are raised to the power of 'gseaParam' before calculation of GSEA enrichment scores.". default = 1. Further comments by fGSEA author at https://github.com/ctlab/fgsea/issues/45
#' @param nPermSimple passed to `fgsea::fgsea()`, from manual: "Number of permutations in the simple fgsea implementation for preliminary estimation of P-values.". default = 50000 in this R package but 1000 by default in fGSEA v1.22.0; we observed much better accuracy in null simulations when increasing this from default 1k to 10k and further minor improvement towards 50k, hence the latter is our default
#' @param gsea_genelist_col optional, only used for `score_type` "custom"
#' @param gsea_scoretype optional, only used for `score_type` "custom"
#' @param random_seed the random seed that is passed to `set.seed()` in order to ensure fgsea results are reproducible. default: 123
#' @return input `genesets` table with results in the "pvalue", "score_type" and "gsea_nes" columns
#' @seealso `test_genesets`
#' @export
test_genesets_gsea = function(genesets, genelist, score_type = NULL, parallel_threads = 1L, gseaParam = 1, nPermSimple = 50000, gsea_genelist_col = NULL, gsea_scoretype = NULL, random_seed = 123) {
  pvalue = effectsize = id_tmp = genes = NULL # fix invisible bindings R package NOTE

  if(!requireNamespace("fgsea", quietly = TRUE)) {
    stop('to perform geneset testing with the GSEA algorithm, the fgsea R package is required (c.f. documentation on the GOAT GitHub page).\nSuggested installation steps;\n1) close RStudio and start it anew (do copy/paste these instructions before closing)\n2) run R command: install.packages("pak")\n3) run R command: pak::pkg_install("ctlab/fgsea", update=FALSE)\n4) test if the installation was successful by loading the R package, run R command: library(fgsea)', call. = FALSE)
  }

  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns source (character type) id (character type), genes (list type)" =
              is.data.frame(genesets) && nrow(genesets) > 0 && all(c("id", "genes") %in% colnames(genesets)) &&
              is.character(genesets$id) && is.list(genesets$genes))
  stopifnot("genelist parameter must be a non-empty data.frame/tibble with column gene (type should match the genes in 'genesets')" =
              length(genelist) > 0 && is.data.frame(genelist) && "gene" %in% colnames(genelist))
  stopifnot("parallel_threads parameter should be a single integer value (0 or higher)" =
              length(parallel_threads) == 1 && is.numeric(parallel_threads) && is.finite(parallel_threads) && parallel_threads >= 0)
  stopifnot("gseaParam parameter should be a single positive numeric value" =
              length(gseaParam) == 1 && is.numeric(gseaParam) && is.finite(gseaParam) && gseaParam > 0)
  stopifnot("nPermSimple parameter should be a single positive numeric value of at least 1000" =
              length(nPermSimple) == 1 && is.numeric(nPermSimple) && is.finite(nPermSimple) && nPermSimple >= 1000)
  stopifnot("score_type parameter must be any of 'pvalue', 'effectsize', 'custom'" = length(score_type) == 1 && score_type %in% c("pvalue", "effectsize", "custom"))
  stopifnot("genelist parameter must have finite non-negative numeric values in the 'pvalue' column when using score_type='pvalue'" =
              (score_type != "pvalue") || ("pvalue" %in% colnames(genelist) && all(is.numeric(genelist$pvalue) & is.finite(genelist$pvalue) & genelist$pvalue >= 0)) )
  stopifnot("genelist parameter must have finite numeric values in the 'effectsize' column when using score_type='effectsize'" =
              (score_type != "effectsize") || ("effectsize" %in% colnames(genelist) && all(is.numeric(genelist$effectsize) & is.finite(genelist$effectsize))) )
  stopifnot("when score_type='custom', genelist parameter must contain a column that matches the column name you specified in parameter gsea_genelist_col (and it should contain only finite numeric values)" =
              (score_type != "custom") || (length(gsea_genelist_col) == 1 && is.character(gsea_genelist_col) && gsea_genelist_col %in% colnames(genelist) &&
                                             all(is.numeric(genelist[,gsea_genelist_col] |> pull()) & is.finite(genelist[,gsea_genelist_col] |> pull())) ) )
  stopifnot("when score_type='custom', gsea_scoretype parameter must be any of 'pos', 'neg', 'std'" =
              (score_type != "custom") || (length(gsea_scoretype) == 1 && gsea_scoretype %in% c('pos', 'neg', 'std')))


  # local helper function
  apply_fgsea = function(gs, gn) {
    pathway = pval = NES = NULL # fix invisible bindings R package NOTE
    # always set.seed prior to this function to enforce the exact same results (on the same computer and input data), regardless
    # of the prior RNG state. Motivating example for this R package where reproducibility of statistical results is essential;
    # Let function A and B both use sample() to draw random numbers.
    # The results from bind_rows(A(), B()) and bind_rows(B(), A()) could be different.
    # Thus, if we do not set.seed() here, every upstream use-case must include a set.seed() call before each function call.
    # Such mistakes are easily made, especially by users to are relatively new to (R) programming, so we enforce it here
    prev_seed <- .Random.seed
    on.exit(.Random.seed <- prev_seed)
    set.seed(random_seed)
    suppressWarnings(fgsea::fgsea(
      pathways = stats::setNames(gs$genes, gs$id), # named geneset list
      stats    = array(gn$score, dimnames = list(gn$gene)), # named gene score array
      scoreType = gsea_scoretype,
      minSize  = 1L, # we already filtered genesets upstream
      maxSize  = 100000L, # we already filtered genesets upstream
      nPermSimple = as.integer(nPermSimple), # Number of permutations in the simple fgsea implementation for preliminary estimation of P-values.  (defaults to 1000 in fGSEA v1.22.0 but here we overwrite this default)
      nproc = as.integer(parallel_threads),
      gseaParam = gseaParam
    )) |>
      tibble::as_tibble() |>
      select(id = pathway, pvalue = pval, gsea_nes = NES)
  }


  # prepare data and parameters for fGSEA
  if(score_type == "pvalue") {
    gsea_scoretype = "pos"
    genelist = genelist |> mutate(score = minlog10_fixzero(pvalue))
  }
  if(score_type == "effectsize") {
    gsea_scoretype = "std"
    genelist = genelist |> mutate(score = effectsize)
  }
  if(score_type == "custom") {
    # use `gsea_scoretype` provided as function parameter
    genelist = genelist |> mutate(score = !!sym(gsea_genelist_col))
  }

  # temp geneset ID that combines source + id for compatibility with input genelists that contain the same ID across 'source'
  genesets$id_tmp = genesets$id
  if("source" %in% colnames(genesets)) {
    genesets$id_tmp = paste(genesets$source, genesets$id)
  }

  result = apply_fgsea(genesets |> select(id = id_tmp, genes), genelist)
  result$score_type = score_type
  if(score_type == "effectsize") {
    result$score_type = c("effectsize_down", "effectsize_up")[1 + (is.finite(result$gsea_nes) & result$gsea_nes > 0)] # if NES > 0, pathway is up-regulated
  }
  # In this function we tested only valid, non-empty genesets that passed filter_genesets() upstream.
  # If the geneset algorithm returned a non-finite pvalue, we consider it "tested but found not signif"
  # (note that leaving NA values as-is would affect subsequent multiple testing correction, which in turn would be a
  # source of bias if many  "no effect whatsoever" genesets are returned as NA pvalues by the algorithm)
  result$pvalue[!is.finite(result$pvalue)] = 1
  result$gsea_nes[!is.finite(result$gsea_nes)] = 0

  genesets |>
    select(-any_of(c("pvalue", "score_type", "gsea_nes"))) |> # remove result columns, if present, prior to left-joining results
    left_join(result, by = c("id_tmp" = "id")) |> # join by temp ID (source + id)
    select(-id_tmp)
}
