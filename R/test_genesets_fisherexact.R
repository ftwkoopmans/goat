
#' Geneset ORA using Fisher-exact test
#'
#' @description In most cases, it's more convenient to call the more generic `test_genesets` function which also deals with multiple-testing correction (per geneset source)
#'
#' It is assumed that the `genesets` and `genelist` parameters are in sync, i.e. `genesets` provided
#' here is the result of the `filter_genesets()` function (using the same `genelist` table)
#'
#' Same as hypergeometric for non-EASE, but slower because stats::fisher.test isn't vectorized
#'
#' Only genesets with at least 1 significant gene are subjected to statistical testing (e.g. NA is returned for genesets without significant genes)
#'
#' @examples \donttest{
#' # note; this example downloads data when first run, and typically takes ~60seconds
#'
#' # store the downloaded files in the following directory. Here, the temporary file
#' # directory is used. Alternatively, consider storing this data in a more permanent location.
#' # e.g. output_dir="~/data/goat" on unix systems or output_dir="C:/data/goat" on Windows
#' output_dir = tempdir()
#'
#' ## first run the default example from test_genesets() to obtain input data
#' datasets = download_goat_manuscript_data(output_dir)
#' genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`
#' genesets_asis = download_genesets_goatrepo(output_dir)
#' genesets_filtered = filter_genesets(genesets_asis, genelist)
#'
#' ## example: same results between Fisher-exact and hypergeometric tests
#' result_hg = test_genesets_hypergeometric(genesets_filtered, genelist, require_nsignif = 3L)
#' result_fe = test_genesets_fisherexact(genesets_filtered, genelist, require_nsignif = 3L)
#' all.equal(result_hg$pvalue, result_fe$pvalue)
#' }
#' @param genesets tibble with genesets, must contain columns 'id', 'ngenes' and 'ngenes_signif'
#' @param genelist tibble with genes, must contain column 'signif'. The number of rows in this table (where signif is not NA)
#' is assumed to be the total number of tested genes, the number of rows where signif==TRUE is assumed the total number of significant genes.
#' @param require_nsignif minimum number of 'signif genes' that overlap with a geneset; `NA` pvalues are returned for genesets with `ngenes_signif <= require_nsignif`.
#' This function 'prefilters' genesets, so beware that this will influence downstream multiple testing correction. Default is 1
#' @param use_ease use a more conservative score coined by DAVID online tools @ https://david.ncifcrf.gov/helps/functional_annotation.html#fisher
#' @return input `genesets` table with results in the "pvalue" column
#' @seealso `test_genesets`
#' @export
test_genesets_fisherexact = function(genesets, genelist, require_nsignif = 1L, use_ease = FALSE) {
  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns ngenes (integer type), ngenes_signif (integer type)" =
              length(genesets) > 0 && is.data.frame(genesets) &&  all(c("ngenes", "ngenes_signif") %in% colnames(genesets)) &&
              is.integer(genesets$ngenes) && is.integer(genesets$ngenes_signif))
  stopifnot("genelist parameter must be a non-empty data.frame/tibble with column 'signif' (logical/boolean type)" =
              length(genelist) > 0 && is.data.frame(genelist) && "signif" %in% colnames(genelist) &&
              is.logical(genelist$signif))
  stopifnot("require_nsignif parameter should be a single positive integer" = length(require_nsignif) == 1 && is.numeric(require_nsignif) && is.finite(require_nsignif) && require_nsignif >= 0)
  stopifnot("use_ease parameter should be a single logical/boolean value" = length(use_ease) == 1 && use_ease %in% c(T,F))

  k_all = sum(!is.na(genelist$signif))     # total number of genes in user's genelist (background)
  k_signif = sum(genelist$signif %in% TRUE)  # total number of significant genes in user's genelist (foreground)
  m = as.integer(use_ease)

  genesets = genesets |> mutate(pvalue = NA_real_)

  # the fisher.test() function in base R is not vectorized, so we have to iterate
  for(i in which(genesets$ngenes_signif >= require_nsignif)) {
    genesets$pvalue[i] = stats::fisher.test(
      matrix(
        c(
          max(0, genesets$ngenes_signif[i] - m),                            # significant genes in geneset G - optional EASE penalty
          genesets$ngenes[i] - genesets$ngenes_signif[i] + m,               # non-signif genes in geneset G ('background' only) + optional EASE penalty
          k_signif - genesets$ngenes_signif[i],                               # significant hits in user input that are not in geneset G
          k_all - k_signif - (genesets$ngenes[i] - genesets$ngenes_signif[i]) # background not in geneset G
        ),
        nrow = 2, byrow = TRUE
      ),
      alternative="greater"
    )$p.value
  }

  return(genesets)
}

