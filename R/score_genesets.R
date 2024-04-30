
#' Compute odds-ratio for each geneset
#'
#' @description
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
#' @param genesets tibble with genesets, must contain columns 'source', 'id', 'ngenes', 'ngenes_signif'
#' @param genelist tibble with genes, must contain columns 'gene', 'signif'
#' @return input `genesets` with results in column "score_oddsratio"
#' @export
score_geneset_oddsratio = function(genesets, genelist) {
  ngenes_signif = ngenes = NULL # fix invisible bindings R package NOTE
  stopifnot(length(genesets) > 0 && is.data.frame(genesets) && all(c("ngenes", "ngenes_signif") %in% colnames(genesets)))
  stopifnot(length(genelist) > 0 && is.data.frame(genelist) && "signif" %in% colnames(genelist))

  k_all = nrow(genelist)                     # total number of genes in user's genelist
  k_signif = sum(genelist$signif %in% TRUE)  # total number of significant genes in user's genelist
  genesets |> mutate(score_oddsratio = (ngenes_signif/ngenes) / (k_signif/k_all))
}



#' Compute a score between -1 and 1 representing the proportion of up- or down-regulated genes for each geneset, weighted by gene effectsizes
#'
#' @description
#' Importantly, the scope/utility of this score is limited to help users post-hoc filter for genesets that contain mostly up/down-regulated genes.
#' However, this might not coincide with the geneset pvalues / significance. For example,
#' genesets may exclusively contain genes with a positive effectsize but at the same time these can all be minor effects/values and thus the geneset
#' is not significant or less significant than other genesets with the exact same "directionality score".
#' For example, genesets may contain both up- and down-regulated genes but still be significant when testing with GOAT and using `score_type='effectsize'`
#'
#' The scores computed with this function can help in post-hoc interpretation of GOAT results to further narrow down all significant genesets to a
#' subset with strong directionality. For example, after `test_genesets()` we can filter the results for
#' A) significant genesets and B) that contain at most N genes and C) that are near-exclusively up/down-regulated.
#' Bringing this all together (also useful for other types of geneset testing, like ORA, score_type="pvalue", etc);
#' `result = test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)`
#' `result = score_geneset_directionality(result, genelist)`
#' `result |> filter(signif & ngenes <= 100 & abs(score_directionality_rank) >= 0.95)`
#'
#' # score definitions;
#' geneset directionality score = weighted mean of respective genes,
#' where gene weights are 1 minus their relative rank in up/down-regulation (depending on negative/positive effectsize)
#' and the value for each gene is -1 or 1 depending on up/down-regulation (sign of effectsize).
#'
#' # pseudocode;
#' 1) gene values and weights
#' A) gene weight between 0 and 1 for the subset of upregulated genes / positive effectsizes;
#' - r = for the subset of genes with effectsize > 0, compute gene rank (1 = highest effectsize, N = smallest effectsize that is greater than zero)
#' - weight = 1 - r/max(r)
#' B) analogous to (A) for the subset of genes with negative effectsize
#' C) result per gene: value = sign of its effectsize, weight = 0 if effectsize 0, otherwise respective weights from (A) or (B)
#'
#' 2) geneset score_directionality = weighted mean over values/weights of respective genes
#'
#' @param genesets tibble with genesets, must contain columns 'source', 'id', 'genes'
#' @param genelist tibble with genes, must contain columns 'gene', 'effectsize'
#' @return input `genesets` with results in 3 columns;
#' `score_directionality_rank` is the weighted gene score where gene values are the sign of their effectsize and weights are linearly proportional to their inverse ranks in the input genelist.
#' `score_directionality_rank2` is similar, but now using rank^2 gene weights to boost the influence of most-important genes in the input genelist.
#' `score_directionality_value` uses the absolute gene effectsizes as gene weights
#' Note that the latter is least robust as it depends on the scaling of input data!
#' @export
score_geneset_directionality = function(genesets, genelist) {
  effectsize = value = weight = source = id = gene = score_directionality = score_directionality_weighted = NULL # fix invisible bindings R package NOTE
  stopifnot(length(genesets) > 0 && is.data.frame(genesets) && all(c("source", "id", "genes") %in% colnames(genesets)) )
  stopifnot(length(genelist) > 0 && is.data.frame(genelist) && all(c("gene", "effectsize") %in% colnames(genelist)) && all(is.numeric(genelist$effectsize) & is.finite(genelist$effectsize)) )

  ### 1) gene scores

  # copy of the genelist
  gl = data.frame(
    gene = genelist$gene,
    effectsize = genelist$effectsize,
    value = sign(genelist$effectsize)
  )

  # compute ranks within subset of up/down
  gl = gl |> arrange(effectsize) # sort by smallest effectsize first (i.e. prio down-regulated)
  gl$rank_down = 1L:nrow(gl)
  gl$rank_down[gl$value != -1] = NA # remove values for not-downregulated

  gl = gl |> arrange(desc(effectsize))
  gl$rank_up = 1L:nrow(gl)
  gl$rank_up[gl$value != 1] = NA

  # compute weights (where 1 = most important gene, near-zero = least important)
  # by rescaling ranks between 0~1 we ensure asymmetric effectsize distributions (proportion of negative/positive values) does not influence relative weights
  gl$rank_down_rescaled = gl$rank_up_rescaled = NA
  if(any(is.finite(gl$rank_down))) { # so we skip max() if there are no non-NA values
    gl$rank_down_rescaled = gl$rank_down / max(gl$rank_down, na.rm = TRUE)
  }
  if(any(is.finite(gl$rank_up))) {
    gl$rank_up_rescaled = gl$rank_up / max(gl$rank_up, na.rm = TRUE)
  }

  gl$weight_down = 1 - gl$rank_down_rescaled
  gl$weight_up = 1 - gl$rank_up_rescaled

  # combine into single weight; i.e. weight_down is only a finite number of effectsize < 0 and vice versa for weight_up
  gl$weight = 0
  gl$weight[is.finite(gl$weight_down)] = gl$weight_down[is.finite(gl$weight_down)]
  gl$weight[is.finite(gl$weight_up)] = gl$weight_up[is.finite(gl$weight_up)]


  ### 2) genesets scores
  # create a new data.frame with only the required data for a minor speedup
  gs = data.frame(source = genesets$source, id = genesets$id, gene = genesets$genes) |>
    tidyr::unchop(cols = gene) |> # convert genesets to long format
    left_join(gl |> select(gene, value, weight, effectsize), by = "gene") |> # add gene-level values/weights
    group_by(source, id) |> # compute geneset scores
    summarise(
      score_directionality_rank = stats::weighted.mean(value, weight),
      score_directionality_rank2 = stats::weighted.mean(value, weight^2),
      score_directionality_value = stats::weighted.mean(value, abs(effectsize)),
      .groups = "drop"
    )
  # scores should always be finite values, but enforce anyway
  gs$score_directionality_rank[!is.finite(gs$score_directionality_rank)] = 0
  gs$score_directionality_rank2[!is.finite(gs$score_directionality_rank2)] = 0
  gs$score_directionality_value[!is.finite(gs$score_directionality_value)] = 0

  # join results into input table
  genesets |>
    select(-any_of(c("score_directionality_rank", "score_directionality_rank2", "score_directionality_value"))) |>
    left_join(gs, by = c("source", "id"))
}
