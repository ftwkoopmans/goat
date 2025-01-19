
#' Test geneset enrichment with the Geneset Ordinal Association Test (GOAT) algorithm
#'
#' @description
#' In most cases, it's more convenient to call the more generic `test_genesets` function which also applies multiple-testing correction (per geneset source) to the geneset p-values computed by this function.
#'
#' This is the canonical geneset test function for GOAT that uses precomputed null distributions that are bundled with the GOAT package
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
#' ### we here compare GOAT with precomputed null distributions against
#' ### a GOAT function that performs bootstrapping to compute null distributions on-demand
#'
#' # apply goat with precomputed null (default) and goat with on-demand bootstrapping
#' result_precomputed = test_genesets(genesets_filtered, genelist, method = "goat",
#'   score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05) |>
#'   # undo sorting by p-value @ test_genesets(), instead sort by stable IDs
#'   arrange(source, id)
#' result_bootstrapped = test_genesets(genesets_filtered, genelist, method = "goat_bootstrap",
#'   score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05, verbose = TRUE) |>
#'   arrange(source, id)
#'
#' # tables should align
#' stopifnot(result_precomputed$id == result_bootstrapped$id)
#' # no missing values
#' stopifnot(is.finite(result_precomputed$pvalue) &
#'           is.finite(is.finite(result_bootstrapped$pvalue)))
#'
#' # compare results
#' plot(result_precomputed$pvalue, result_bootstrapped$pvalue)
#' abline(0, 1, col=2)
#'
#' plot(minlog10_fixzero(result_precomputed$pvalue),
#'      minlog10_fixzero(result_bootstrapped$pvalue))
#' abline(0, 1, col=2)
#'
#' summary(minlog10_fixzero(result_precomputed$pvalue) -
#'         minlog10_fixzero(result_bootstrapped$pvalue))
#' }
#' @param genesets genesets data.frame, must contain columns; "source", "id", "genes", "ngenes"
#' @param genelist genelist data.frame, must contain columns "gene" and "pvalue"/"effectsize" (depending on parameter `score_type`)
#' @param score_type how to compute gene scores?
#' Option "pvalue" uses values from the pvalue column in `genelist` in a one-way test for enrichment; lower p-value is better
#' Option "effectsize" uses values from the effectsize column in `genelist` in a two-way test for enrichment; is a geneset enriched in either down- or up-regulated genes?
#' Option "effectsize_abs" uses values from the effectsize column in `genelist` in a one-way test for enrichment; is a geneset enriched when testing absolute effectsizes?
#' Option "effectsize_up" uses values from the effectsize column in `genelist` in a one-way test for enrichment; is a geneset enriched in up-regulated genes?  (i.e. positive effectsize)
#' Option "effectsize_down" uses values from the effectsize column in `genelist` in a one-way test for enrichment; is a geneset enriched in down-regulated genes?  (i.e. negative effectsize)
#' @return input `genesets` table with results in the "pvalue", "score_type" columns.
#' "zscore" column:
#' A standardized z-score is computed from geneset p-values + effectsize direction (up/down) if tested.
#' Importantly, we here return standardized z-scores because the GOAT geneset score (mean of gene scores) is relative to the respective geneset-size-matched null distributions (a skewed normal)!
#' In contrast, the standardized z-scores are comparable between genesets (as are the pvalues).
#'
#' Only if either (or both) the effectsize-up/down was tested, the direction of regulation has been tested (effectsize_abs and pvalue score types are agnostic to up/down regulation).
#' So when score_type was set to any of effectsize/effectsize_down/effectsize_up, the z-scores are negative values in case the "score_type" output column is "effectsize_down".
#' @seealso `test_genesets`
#' @export
test_genesets_goat_precomputed = function(genesets, genelist, score_type) {
  N = NULL
  validate_goat_scoretype(score_type)
  validate_goat_genelist(genelist, score_type)
  validate_goat_genesets(genesets, nrow(genelist))

  goat_null = goat::goat_nulldistributions

  ### compute null distribution parameters
  geneset_usize = sort(unique(genesets$ngenes), decreasing = FALSE)
  goat_nulldistributions_subset = goat_null |> filter(N == nrow(genelist))
  # must have exactly 1 matching row, otherwise there is no precomputed null for this genelist size
  stopifnot("genelist length not available in precomputed null distributions" = nrow(goat_nulldistributions_subset) == 1)
  # double-check output
  nulldistribution_parameters = goat_apply_precomputed_null(usize = geneset_usize, par = goat_nulldistributions_subset)

  # perform geneset testing using the null distribution parameters obtained above
  test_genesets_goat(genesets, genelist, score_type, nulldistribution_parameters)
}



#' Variant of the main GOAT function `test_genesets_goat_precomputed` that does not use previously prepared parameters
#'
#' In typical use-cases, one applies `test_genesets()` instead with parameter `method="goat"` ,
#' which in turn will use `test_genesets_goat_precomputed` (and not this function).
#'
#' Optionally, use this function to ignore precomputed/bundled null distribution estimates and perform new permutations and function fitting
#' (e.g. because you want to test the effect of huge niter parameter, but beware of RAM requirements)
#'
#' @param genesets see `test_genesets_goat_precomputed`
#' @param genelist see `test_genesets_goat_precomputed`
#' @param score_type see `test_genesets_goat_precomputed`
#' @param niter integer number of bootstrap iterations; at least 10000, at most 5000000
#' @param verbose boolean, create debug plots
#' @return see `test_genesets_goat_precomputed`
#' @export
test_genesets_goat_fitfunction = function(genesets, genelist, score_type, niter = 500000, verbose = FALSE) {
  validate_goat_scoretype(score_type)
  validate_goat_genelist(genelist, score_type)
  validate_goat_genesets(genesets, nrow(genelist))
  validate_goat_niter(niter)
  validate_goat_verbose(verbose)

  geneset_usize = sort(unique(genesets$ngenes))
  l = goat_nulldistribution_function(nrow(genelist), niter = niter, return_fit_objects = FALSE, verbose = verbose)
  nulldistribution_parameters = goat_apply_precomputed_null(usize = geneset_usize, par = l)

  # perform geneset testing using the null distribution parameters obtained above
  test_genesets_goat(genesets, genelist, score_type, nulldistribution_parameters)
}



#' Naive GOAT variant where we estimate null parameters for each geneset size independently
#'
#' In typical use-cases, one applies `test_genesets()` instead with parameter `method="goat"` ,
#' which in turn will use `test_genesets_goat_precomputed` (and not this function).
#'
#' @param genesets see `test_genesets_goat_precomputed`
#' @param genelist see `test_genesets_goat_precomputed`
#' @param score_type see `test_genesets_goat_precomputed`
#' @param niter integer number of bootstrap iterations; at least 10000, at most 5000000
#' @param verbose boolean, create debug plots
#' @return see `test_genesets_goat_precomputed`
#' @export
test_genesets_goat_bootstrap = function(genesets, genelist, score_type, niter = 500000, verbose = FALSE) {
  validate_goat_scoretype(score_type)
  validate_goat_genelist(genelist, score_type)
  validate_goat_genesets(genesets, nrow(genelist))
  validate_goat_niter(niter)
  validate_goat_verbose(verbose)

  # for each unique geneset size, use bootstrapping to generate a null distribution then fit a skew-normal
  # result is data.frame with estimated parameters (of each fitted skew-normal)
  nulldistribution_parameters = goat_nulldistribution_independent(
    geneset_sizes = genesets$ngenes,
    # the order does not matter; this is the overall distribution of gene scores to sample from when computing geneset scores
    gene_scores = rankscore_fixed_order(nrow(genelist)),
    niter = niter,
    verbose = verbose
  )

  # perform geneset testing using the null distribution parameters obtained above
  test_genesets_goat(genesets, genelist, score_type, nulldistribution_parameters)
}



#' Compute the parameters of a skew-normal distribution, given a previously fitted polynomial function over geneset sizes
#'
#' @param usize integer vector of geneset sizes
#' @param par result from `goat_nulldistribution_function`
#' @return data.frame with columns; size, mu, sigma, xi
#' @noRd
goat_apply_precomputed_null = function(usize, par) {
  stopifnot(is.list(par) && length(par$sigma_0) == 1)
  f_xi = function(x,l) l$xi_0 + l$xi_1 * x + l$xi_2 * x^2 + l$xi_3 * x^3 + l$xi_4 * x^4 + l$xi_5 * x^5 + l$xi_6 * x^6 + l$xi_7 * x^7 + l$xi_8 * x^8 + l$xi_9 * x^9 + l$xi_10 * x^10 + l$xi_11 * x^11
  f_sigma = function(x,l) l$sigma_0 + l$sigma_1 * x + l$sigma_2 * x^2 + l$sigma_3 * x^3 + l$sigma_4 * x^4 + l$sigma_5 * x^5 + l$sigma_6 * x^6 + l$sigma_7 * x^7 + l$sigma_8 * x^8 + l$sigma_9 * x^9

  data.frame(
    size = usize,
    mu = rep(par$mu, length(usize)),
    sigma = f_sigma(log(usize), par),
    xi = pmax(1, exp(f_xi(log(usize), par)) )
  )
}



#' From effectsize and pvalue in the genelist table to a matrix of 'rank transformed gene scores'
#'
#' @description
#' Depending on the score_type, we sort the genelist by pvalue or effectsize (up/down/absolute) and
#' then compute gene scores. Note that it is common for genelists to have many genes with `pvalue=1`,
#' which result in ties when sorting. We break ties for pvalue sorting by absolute effectsizes
#' (if available) and we break ties for (up/down/absolute) effectsizes by pvalues (if available).
#' Finally, we break any remaining ties by sorting the gene column to ensure stable sorting behaviour
#' under all circumstances.
#'
#' Technical note; after upstream input validation, we're assured this will yield a non-empty matrix with at least 1 column
#'
#' @param genelist see `test_genesets_goat_precomputed`
#' @param score_type see `test_genesets_goat_precomputed`
#' @return gene*score_type matrix where rows are aligned with the rows of the genelist table
#' @noRd
goat_testgene_score_matrix = function(genelist, score_type) {
  has_pvalue = "pvalue" %in% colnames(genelist)
  has_effectsize = "effectsize" %in% colnames(genelist)
  genelist_scores = NULL

  # effectsize_up
  if(any(score_type %in% c("effectsize", "effectsize_up"))) {
    if(has_pvalue) {
      x = rankscore(genelist, genelist$effectsize, -1 * genelist$pvalue, genelist$gene, colname = "score")
    } else {
      x = rankscore(genelist, genelist$effectsize, genelist$gene, NULL, colname = "score")
    }
    genelist_scores = cbind(genelist_scores, effectsize_up = x$score)
  }

  # effectsize_down; analogous to above, but with -1 * effectsize to rank "lower is better"
  if(any(score_type %in% c("effectsize", "effectsize_down"))) {
    if(has_pvalue) {
      x = rankscore(genelist, -1 * genelist$effectsize, -1 * genelist$pvalue, genelist$gene, colname = "score")
    } else {
      x = rankscore(genelist, -1 * genelist$effectsize, genelist$gene, NULL, colname = "score")
    }
    genelist_scores = cbind(genelist_scores, effectsize_down = x$score)
  }

  # effectsize_abs; analogous, but with absolute effectsize
  if(any(score_type %in% c("effectsize_abs"))) {
    if(has_pvalue) {
      x = rankscore(genelist, abs(genelist$effectsize), -1 * genelist$pvalue, genelist$gene, colname = "score")
    } else {
      x = rankscore(genelist, abs(genelist$effectsize), genelist$gene, NULL, colname = "score")
    }
    genelist_scores = cbind(genelist_scores, effectsize_abs = x$score)
  }

  # pvalue
  if(any(score_type %in% c("pvalue"))) {
    if(has_effectsize) {
      x = rankscore(genelist, -1 * genelist$pvalue, abs(genelist$effectsize), genelist$gene, colname = "score")
    } else {
      x = rankscore(genelist, -1 * genelist$pvalue, genelist$gene, NULL, colname = "score")
    }
    genelist_scores = cbind(genelist_scores, pvalue = x$score)
  }

  stopifnot("unknown 'score_type' parameter" = !is.null(genelist_scores))
  return(genelist_scores)
}



#' Implementation of the GOAT algorithm
#'
#' @param genesets see `test_genesets_goat_precomputed`
#' @param genelist see `test_genesets_goat_precomputed`
#' @param score_type see `test_genesets_goat_precomputed`
#' @param nulldistribution_parameters result from goat_nulldistribution_independent
#' @noRd
test_genesets_goat = function(genesets, genelist, score_type, nulldistribution_parameters) {
  stopifnot("precomputed null distributions contain non-finite values" = is.finite(nulldistribution_parameters$size) &
              is.finite(nulldistribution_parameters$mu) &
              is.finite(nulldistribution_parameters$sigma) &
              is.finite(nulldistribution_parameters$xi))

  ### 1) from genelist table + score_type to gene scores
  # matrix where columns are 'score types' (e.g. "effectsize_up" and "effectsize_down") and rows are genes (aligned with genelist table)
  genelist_scores = goat_testgene_score_matrix(genelist, score_type)

  ### 2) compute geneset scores using C++ helper for fast aggregation
  ul_gs_index = rep(seq_len(nrow(genesets)), genesets$ngenes)
  ul_gs_geneindex = match(unlist(genesets$genes, use.names = FALSE, recursive = FALSE), genelist$gene)
  # double-check that this is the filtered geneset input; all genes in the genesets table must match the genelist table
  stopifnot("all genes in the genesets table must be available in the genelist table  (did you forget to apply the filter_genesets() function?)" = !anyNA(ul_gs_geneindex))
  # importantly, rcpp_gene_to_geneset_scores returns sum scores not the mean !
  gs_scores = rcpp_gene_to_geneset_scores(n_gs = nrow(genesets), gs_index = ul_gs_index, gs_geneindex = ul_gs_geneindex, gene_score = genelist_scores)
  # note that we still have to divide the geneset scores computed in C++ by the geneset size (there is only summation in this C++ function)
  genesets$score_type = colnames(genelist_scores)[1]
  if(ncol(genelist_scores) == 1) {
    genesets$score = gs_scores[,1] / genesets$ngenes # !! still have to compute mean
  } else {
    stopifnot(ncol(genelist_scores) == 2) # currently, there are at most 2 geneset score types computed at the same time so we can use this efficient code
    genesets$score = pmax(gs_scores[,1], gs_scores[,2]) / genesets$ngenes # !! still have to compute mean
    genesets$score_type = colnames(genelist_scores)[1 + (gs_scores[,2] > gs_scores[,1])] # use colname matching highest score
  }

  ### 3) geneset score to pvalue
  # find the respective null distribution parameters (mu,sigma,xi @ skew-normal) by matching by geneset size (gene count)
  # matching never yields NA because we already checked all unique 'ngenes' are present in 'nulldistribution_parameters' upstream
  i = match(genesets$ngenes, nulldistribution_parameters$size)
  genesets$pvalue = psnorm_upper_tail(
    genesets$score,
    mean = nulldistribution_parameters$mu[i],
    sd = nulldistribution_parameters$sigma[i],
    xi = nulldistribution_parameters$xi[i]
  )
  genesets$score = NULL # no longer need the score column, remove it

  # finally, adjust for multiple scores (if applicable)
  if(ncol(genelist_scores) > 1) {
    genesets$pvalue = pmin(1, genesets$pvalue * ncol(genelist_scores))
  }

  # A standardized z-score is computed from geneset p-values + effectsize direction (up/down) if tested.
  genesets$zscore = stats::qnorm(genesets$pvalue / 2, lower.tail = FALSE)
  if(any(score_type %in% c("effectsize", "effectsize_down", "effectsize_up"))) {
    genesets$zscore = genesets$zscore * c(-1,1)[1L + as.integer(genesets$score_type == "effectsize_up")]
  }

  return(genesets)
}



#' helper function to validate input parameters for GOAT computation functions
#' @param x object to validate
#' @noRd
validate_goat_scoretype = function(x) {
  stopifnot("score_type parameter must be any of 'pvalue', 'effectsize', 'effectsize_abs', 'effectsize_up', 'effectsize_down'" = length(x) == 1 && x %in% c("pvalue", "effectsize", "effectsize_up", "effectsize_down", "effectsize_abs"))
}



#' helper function to validate input parameters for GOAT computation functions
#' @param x object to validate
#' @noRd
validate_goat_niter = function(x) {
  stopifnot("niter parameter must be a single integer between 10000 and 5000000" = length(x) == 1 && is.numeric(x) && is.finite(x) && x >= 10000L && x <= 5000000)
}



#' helper function to validate input parameters for GOAT computation functions
#' @param x object to validate
#' @noRd
validate_goat_verbose = function(x) {
  stopifnot("verbose parameter must be a single boolean (TRUE or FALSE)" = length(x) == 1 && x %in% c(TRUE, FALSE))
}



#' helper function to validate input parameters for GOAT computation functions
#' @param x object to validate
#' @param score_type score type
#' @noRd
validate_goat_genelist = function(x, score_type) {
  stopifnot("genelist parameter must be a non-empty data.frame/tibble and contain a 'gene' column" =
              length(x) > 0 && is.data.frame(x) && "gene" %in% colnames(x) )
  stopifnot("genelist parameter must be a table of at least 100 genes" = nrow(x) >= 100)
  stopifnot("genelist parameter must be a table of at most 50000 genes" = nrow(x) <= 50000)
  stopifnot("score_type parameter must be a single string" = length(score_type) == 1 && is.character(score_type))
  # depending on score_type, check for valid data in pvalue and effectsize columns
  stopifnot("genelist parameter must have finite non-negative numeric values in the 'pvalue' column when using score_type='pvalue'" =
              (score_type != "pvalue") || ("pvalue" %in% colnames(x) && all(is.numeric(x$pvalue) & is.finite(x$pvalue) & x$pvalue >= 0)) )
  stopifnot("genelist parameter must have finite numeric values in the 'effectsize' column when using a 'score_type' based on effectsizes" =
              (!score_type %in% c("effectsize", "effectsize_abs", "effectsize_up", "effectsize_down")) || ("effectsize" %in% colnames(x) && all(is.numeric(x$effectsize) & is.finite(x$effectsize))) )
}



#' helper function to validate input parameters for GOAT computation functions
#' @param x object to validate
#' @param genelist_N genelist length
#' @noRd
validate_goat_genesets = function(x, genelist_N) {
  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns source, id, genes, ngenes" =
              length(x) > 0 && is.data.frame(x) && nrow(x) > 0 && all(c("source", "id", "genes", "ngenes") %in% colnames(x)) )
  stopifnot("genesets parameter must contain integer type column 'ngenes' (these are used for indexing, so convert using as.integer upstream if needed)" =
              is.integer(x$ngenes) )
  ngene_range = range(x$ngenes)
  stopifnot("the maximum geneset size ('ngenes' column in genesets table) cannot be larger than half the genelist length (rows in genelist table).\nDid you forget to apply the filter_genesets() function?" = ngene_range[2] <= genelist_N)
  stopifnot("genesets parameter contains genesets with a size (gene count) that is not compatible with GOAT's null distributions. Supported geneset sizes are at least 10 and at most half the genelist length.\nTo fix this, configure the function filter_genesets() with parameters min_overlap=10 and max_overlap_fraction=0.5\n(additionally, the recommended maximum absolute geneset size for typical analyses of the GO database is 1500~2000 , which you can configure with parameter max_overlap=1500 in the filter_genesets function)" =
              ngene_range[1] >= 10 && ngene_range[2] <= ceiling(genelist_N/2)) # we only precomputed / fit to up 50% of genelist size
}
