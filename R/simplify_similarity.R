
#' Reduce the set of significant genesets to a minimum
#'
#' Analyses are performed independently per 'source' of genesets. The result of this function is the geneset table with a newly appended column 'signif_and_reduced'
#'
#' @param clusters results from `cluster_genesets()`
#' @param simscore_threshold similarity score (0~1) that is required to consider one geneset to be a "parent term" of another.
#' Setting a lower value will yield fewer genesets / stronger summarization.
#' Typical settings for this parameter are 0.8~0.99 (0.9 is default)
#' @param universe_fraction discard genesets that cover more than X fraction of all genes in the universe (unique set of genes covered by all significant genesets).
#' Setting this to 0.25 will deprioritize genesets that cover 25% of all genes (in significant genesets).
#' This prevents very generic GO terms like "protein-containing complex" to be included in results.
#' Typical settings for this parameter are 0.1~0.5  (0.25 is default)
#' @param signifgenes_fraction the minimum fraction of "foreground genes" ('genes_signif' column) found across all significant genesets that should be covered by the reduced geneset collection.
#' This parameter doesn't do anything if there are fewer than 5 "foreground genes" alltogether.
#' Typical settings for this parameter are 0.75~0.95 (0.9 is default)
#' @export
reduce_genesets = function(clusters, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9) {
  ngenes_signif = genes_signif = rescue_count = ngenes = rescue_score = ngenes_input = NULL # fix invisible bindings R package NOTE
  # init result column
  clusters$genesets = clusters$genesets |> mutate(signif_and_reduced = signif)

  for(src in names(clusters$similarity)) { # src = "GO_CC"
    src_rows = clusters$genesets$signif %in% TRUE & clusters$genesets$source == src
    result_subset = clusters$genesets |> filter(src_rows)
    if(nrow(result_subset) < 2) next # nothing to do
    x = clusters$similarity[[src]] # geneset*geneset similarity matrix where scores 0~1 indicate if col j is a parent of row i
    diag(x) = 0 # remove self-similarities
    universe = n_distinct(unlist(result_subset$genes)) # count unique genes covered by all significant genesets

    # goal; ignore terms that are parents of more than X% of the universe
    ids_ignore = NULL
    # iterate columns; values on rows indicate if j is a parent of <row>
    for(j in 1:ncol(x)) {
      # find all children of j, but disregard those that are already on the ignore list.
      # this prevents scenarios where 3 highly overlapping terms, that cover a large part of the universe,
      # are all removed under the assumption that one of the other terms is a child term that'll be retained
      j_children = x[ ! rownames(x) %in% ids_ignore, j] >= simscore_threshold
      if(sum(j_children) > 0) {
        # from the geneset table obtained through test_genesets(), collect the list of genes for genesets that are in geneset j or any of its child terms
        g = result_subset$genes[ result_subset$id %in% colnames(x)[c(j, which(j_children))] ]
        n = n_distinct(unlist(g)) # count unique genes
        if(n / universe > universe_fraction) {
          ids_ignore = c(ids_ignore, colnames(x)[j])
        }
      }
    }

    # subset the similarity matrix
    X = x[ ! rownames(x) %in% ids_ignore, ! colnames(x) %in% ids_ignore, drop = FALSE]
    # find genesets that do not have parents -->> we retain these as representatives of all genesets
    X = X > simscore_threshold
    ids_retain = rownames(X)[rowSums(X) == 0]

    # deal with edge-cases where after pruning terms, less than X% of the universe is covered.
    # here recursively reinstate ids_ignore that cover the most 'missing genes'
    foreground_universe = unique(unlist(result_subset$genes_signif))
    if(length(foreground_universe) >= 5) {
      foreground_covered = unique(unlist(result_subset$genes_signif[result_subset$id %in% ids_retain]))
      while(length(foreground_covered) / length(foreground_universe) < signifgenes_fraction) {
        foreground_universe_missing = setdiff(foreground_universe, foreground_covered)
        # order the ignore IDs by the fraction of genelist they cover (i.e. include terms in order of 'purity')
        tmp = result_subset |> filter(id %in% ids_ignore & ngenes_signif >= 3) |> mutate(rescue_count = lengths(lapply(genes_signif, intersect, foreground_universe_missing)), rescue_score = rescue_count / ngenes) |> arrange(desc(rescue_score), ngenes_input)
        # # order the ignore IDs by the absolute number of significant genes they'll contribute
        # tmp = result_subset |> filter(id %in% ids_ignore) |> mutate(rescue_score = lengths(lapply(genes_signif, intersect, foreground_universe_missing))) |> arrange(desc(rescue_score), ngenes_input)
        ids_retain = c(ids_retain, tmp$id[1])
        foreground_covered = unique(unlist(result_subset$genes_signif[result_subset$id %in% ids_retain]))
      }
    }

    # finally, store results by updating clusters$genesets
    clusters$genesets$signif_and_reduced[src_rows] = clusters$genesets$id[src_rows] %in% ids_retain
  }

  return(clusters$genesets)
}



#' cluster significant genesets from `test_genesets()` by geneset similarity (separately for each 'geneset source')
#'
#' @param x results from `test_genesets()`
#' @param genelist should be the same as provided to `test_genesets()`
#' @param hclust_method hierarchical clustering method, any of; 'ward.D', 'ward.D2' (default), 'single', 'complete', 'average'
#' @export
cluster_genesets = function(x, genelist, hclust_method = "ward.D2") {
  x = validate_genesets(x)
  genelist = validate_genelist(genelist)
  # in addition to standard geneset validation, also check that this is a 'genelist' result; we need to know which are significant
  stopifnot("geneset table should contain a 'signif' column with boolean values (i.e. the result from test_genesets() )" = (
    "signif" %in% colnames(x) && is.logical(x$signif) # allow NA values @ signif
  ))
  stopifnot("'hclust_method' parameter must be any of; 'ward.D', 'ward.D2', 'single', 'complete', 'average'" =
              length(hclust_method) == 1 && hclust_method %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average'))


  score_type = unique(x$score_type)
  # default to a distance matrix based on pvalue-based gene scores
  if(length(score_type) == 0) {
    score_type = "pvalue"
  }

  stopifnot("score_type column in input table can only contain these values; 'effectsize', 'effectsize_up', 'effectsize_down', 'effectsize_abs', 'effectsize_3way', 'pvalue', 'pvalue_effectsize'" = (
    length(score_type) > 0 && all(score_type %in% c("effectsize", "effectsize_up", "effectsize_down", "effectsize_abs", "effectsize_3way", "pvalue", "pvalue_effectsize"))
  ))

  ### genelist scores, analogous to test_genesets_goat()
  # matrix where columns are 'score types' (e.g. "effectsize_up" and "effectsize_down") and rows are genes (aligned with genelist table)
  genelist_scores = goat_testgene_score_matrix(genelist, score_type)

  # init result variables
  rows = x$signif %in% TRUE
  similarity = hc_row = hc_col = list()

  for(src in sort(unique(x$source))) { # src = "GO_CC"
    rows = x$source == src & x$signif %in% TRUE # use %in% to guard against NA values / 'not tested'
    N = sum(rows)
    # nothing to do, go next
    if(N == 0) next
    # don't compute distance matrix nor clustering if there is only 1 geneset
    if(N == 1) {
      similarity[[src]] = matrix(1, nrow = 1, ncol = 1, dimnames = list(x$id[rows], x$id[rows]))
      next
    }

    # compute distance matrix
    y = x |> filter(rows)
    m = genesets_as_matrix(y)
    sim = geneset_similarity_matrix(mat = m$mat, mat_row_geneid = m$rowid, mat_col_gsid = m$colid, genelist_geneid = genelist$gene, genelist_scores = genelist_scores, weighted = TRUE)
    # if(symmetric_similarity_matrix) {
    #   sim = (t(sim) + sim) / 2
    # }
    similarity[[src]] = sim

    hc_row[[src]] = stats::hclust(stats::dist(sim), method = hclust_method)
    hc_col[[src]] = stats::hclust(stats::dist(t(sim)), method = hclust_method)
  }

  # settings as string
  settings = sprintf("simplify_genesets(hclust_method='%s')", hclust_method)
  attr(x, "settings") <- c(attr(x, "settings"), settings)

  return(list(genesets = x, similarity = similarity, hc_row = hc_row, hc_col = hc_col))
}



#' utility function that converts a genesets table into a sparse gene*geneset identity matrix
#'
#' @param genesets genesets tibble. Must contain columns; "id", "genes", "ngenes"
genesets_as_matrix = function(genesets) {
  stopifnot(c("id", "genes", "ngenes") %in% colnames(genesets))
  # unlist the geneset table (row = geneset, genes are packed in a list column)
  ul_id = rep(genesets$id, genesets$ngenes)
  ul_gene = unlist(genesets$genes, recursive = F, use.names = F)
  # ul_score = genelist$score_pval[match(ul_gene, genelist$gene)]  # this strategy is much slower on huge datasets
  # tmp = aggregate(ul_score, by = list(ul_id), FUN=sum)             # this strategy is much slower on huge datasets
  ugene = unique(ul_gene)
  # return a sparse gene*geneset identity matrix
  return( list(mat = Matrix::sparseMatrix(i = match(ul_gene, ugene),
                                          j = match(ul_id, genesets$id)),
               rowid = ugene,
               colid = genesets$id) )
}



#' geneset similarity matrix computation
#'
#' importantly, if multiple genelist_scores are provided, we use the gene score type for geneset j that yields the highest sum score
#' @param mat numeric matrix
#' @param mat_row_geneid for mat rows; gene identifiers
#' @param mat_col_gsid for mat columns; geneset identifiers
#' @param genelist_geneid same identifiers as `mat_row_geneid`
#' @param genelist_scores a matrix where columns are various gene score types
#' @param weighted boolean value
geneset_similarity_matrix = function(mat, mat_row_geneid, mat_col_gsid, genelist_geneid, genelist_scores, weighted = TRUE) {

  ### 1) find the score for each row in the gene*geneset matrix via lookup in the genelist table
  # geneset score; gene score vector along all rows * identity matrix = matrix with gene score values in respective genesets -->> sum of gene scores / ngenes = geneset scores
  i = match(mat_row_geneid, genelist_geneid)
  stopifnot("all genes in the genelist table must be present in the geneset matrix" = !anyNA(i))

  ### 2) asymmetric distance matrix;
  # gene identity vector for geneset g * identity matrix = matrix with 'gene on row x is present in both genesets'
  Ngs = ncol(mat)
  d = matrix(0, nrow = Ngs, ncol = Ngs, dimnames = list(mat_col_gsid, mat_col_gsid))
  j_scoretypes = rep(0, ncol(genelist_scores))
  # iterate columns in the genes*genesets identity matrix
  for(j in seq_len(Ngs)) {
    if(weighted) {
      ## A) for current geneset j, compute the 'gene score type' (e.g. rank transformed p-values or effectsizes)
      ##    that yields the highest geneset score  =  score type most relevant for this geneset -->> use that for similarity computation downstream
      for(index_scoretype in 1L:ncol(genelist_scores)) {
        j_scoretypes[index_scoretype] = sum(mat[,j] * genelist_scores[i, index_scoretype])
      }
      index_scoretype = which.max(j_scoretypes)

      ## B) weighted overlap between geneset j (current loop iteration) and all genesets (columns in 'mat')
      j_score = mat[,j] * genelist_scores[i, index_scoretype]    # genes without a score in geneset j are set to zero score through multiplication with identity matrix 'mat'
      j_overlap_score = Matrix::colSums(mat * j_score)            # sum of gene scores in overlapping set of genes present in geneset j and all other genesets
      d[j,] = j_overlap_score / j_overlap_score[j]                # divide by the total sum of gene scores in geneset j; result = proportion of gene scores in geneset j identified in other genesets
    } else {
      # not weighted = simple case; count number of overlapping genes
      j_overlap_score = Matrix::colSums(mat * mat[,j])
      d[j,] = j_overlap_score / j_overlap_score[j] # overlap count for geneset_ij, divided by respective geneset size = proportion of overlap
    }
  }

  return(d)
}



#' plot the geneset similarity matrix as a heatmap
#'
#' @examples \dontrun{
#'   # try various color palettes
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Viridis", rev = F))
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Inferno", rev = F))
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Lajolla", rev = T))
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Mako", rev = F))
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Turku", rev = T))
#'   plot_heatmap(plot_heatmapclusters, heatmap_colors = hcl.colors(100, "Grays", rev = T))
#' }
#' @param x result from `cluster_genesets()`
#' @param output_dir full path to output directory. Set to NA to directly show the figures instead of writing them to file
#' @param heatmap_colors a vector of 100 colors to be used for the heatmap (101 breaks are computed between 0 and the max value in the distance matrix)
#' @param fontsize parameter sent to pheatmap::pheatmap(); control the size of labels in the plot, defaults to 10. Note that you can also change the plot device size, see examples
#' @export
plot_heatmap = function(x, output_dir = NA, heatmap_colors = grDevices::hcl.colors(100, "Viridis", rev = F), fontsize = 10) {
  check_dependency("pheatmap", "plotting heatmaps")
  stopifnot("parameter output_dir must be a full path for the output files, this directory must already exists" =
              length(output_dir) == 1 && (is.na(output_dir) || (is.character(output_dir) && nchar(output_dir) > 4 && dir.exists(output_dir)) ) )
  stopifnot("heatmap_colors parameter must be a vector of 100 html color-codes" = length(heatmap_colors) == 100 && is.character(heatmap_colors) && all(grepl("^#[0-9a-zA-Z]{6}$", heatmap_colors)))
  stopifnot("fontsize parameter must be a positive number" = length(fontsize) == 1 && is.numeric(fontsize) && is.finite(fontsize) && fontsize > 0)
  if(length(x) == 0 || !is.list(x) || !"similarity" %in% names(x) || length(names(x$similarity)) == 0) {
    cat("empty data provided to heatmap plot function\n")
    return()
  }
  stopifnot("parameter x must be a result from the cluster_genesets() function; cannot find a geneset table" =
              "genesets" %in% names(x) && is.data.frame(x$genesets) && all(c("id", "source", "name", "ngenes") %in% colnames(x$genesets)))

  fontsize = as.integer(round(fontsize))

  score_types = NULL
  if("score_type" %in% colnames(x$genesets)) {
    score_types = sort(stats::na.omit(unique(x$genesets$score_type)))
  }

  for(src in names(x$similarity)) {
    hm_data = x$similarity[[src]]
    if(length(hm_data) == 0 || !is.matrix(hm_data)) {
      next
    }
    if(!all(is.finite(hm_data))) {
      stop(paste("invalid similarity object; distance matrix contains non-finite values @", src))
    }
    mat_names = rownames(hm_data)
    if(length(mat_names) != nrow(hm_data)) {
      stop(paste("invalid similarity object; missing rownames in the distance matrix @", src))
    }
    i = match(mat_names, x$genesets$id)
    if(anyNA(i)) {
      stop(paste("invalid similarity object; some rownames in the distance matrix don't match the geneset table @", src))
    }
    mat_ngenes = x$genesets$ngenes[i]
    mat_labels = string_trunc_right(x$genesets$name[i], width = 50)
    mat_names_scoretype = rep("", length(i))

    if(length(score_types) > 0) {
      mat_names_scoretype = x$genesets$score_type[i]
    }


    hc_row = hc_col = NULL
    if("hc_row" %in% names(x) && is.list(x$hc_row) && src %in% names(x$hc_row)) {
      hc_row = x$hc_row[[src]]
    }
    if("hc_col" %in% names(x) && is.list(x$hc_col) && src %in% names(x$hc_col)) {
      hc_col = x$hc_col[[src]]
    }

    is_symmetric = isSymmetric(hm_data)
    do_cluster = nrow(hm_data) > 2 && length(hc_row) > 0
    if(!do_cluster) {
      hc_row = hc_col = FALSE
    }


    ## manually bin the geneset size data into 10 distinct values and assign colors (i.e. don't use pheatmap defaults)
    # breaks from 0 to max value
    # nbreaks = 11
    # mat_ngenes_breaks = seq(from = min(c(10, mat_ngenes)), to = max(c(mat_ngenes, 50)), length.out = nbreaks)
    mat_ngenes_breaks = pretty(c(min(c(10, mat_ngenes)), max(c(mat_ngenes, 50))), n = 8)
    nbreaks = length(mat_ngenes_breaks)
    mat_ngenes_colors = stats::setNames(
      grDevices::hcl.colors(nbreaks - 1, palette = "Grays", rev = T),
      mat_ngenes_breaks[-1]
    )
    # bin / color index (goes to nbreaks-1 max)
    mat_ngenes_binned = rep(1L, length(mat_ngenes))
    for(i in 1:(nbreaks-1)) {
      mat_ngenes_binned[mat_ngenes > mat_ngenes_breaks[i]] = i
    }

    annot_row = data.frame(
      `gene count` = mat_ngenes_breaks[mat_ngenes_binned],
      row.names = rownames(hm_data),
      check.names = F
    )
    annot_clr = list(
      `gene count` = mat_ngenes_colors
    )
    if(length(score_types) > 0 && any(mat_names_scoretype != "")) {
      annot_row$score_type = as.character(mat_names_scoretype)
      annot_clr$score_type = stats::setNames(grDevices::rainbow(length(score_types)), score_types)
    }

    heatmap_breaks = seq(from = 0, to = 1, length.out = 101)

    f = NA
    if(!is.na(output_dir)) {
      f = sprintf("%s/goat_heatmap_%s.pdf", output_dir, gsub("[^a-zA-Z0-9 _-]+", "", src))
      if(file.exists(f)) {
        file.remove(f)
        if(file.exists(f)) {
          stop(paste("failed to overwrite existing file when creating a heatmap plot:", f))
        }
      }
    }

    p = pheatmap::pheatmap(
      hm_data,
      cluster_rows = hc_row,
      cluster_cols = hc_col,
      labels_row = mat_labels,
      labels_col = mat_labels,
      fontsize = fontsize - 3L,
      fontsize_col = fontsize,
      fontsize_row = fontsize,
      annotation_row = annot_row,
      annotation_colors = annot_clr,
      annotation_legend = TRUE,
      breaks = heatmap_breaks,
      color = heatmap_colors,
      cellwidth = 10,
      cellheight = 10,
      filename = f,
      main = sprintf(
        "%s - %d genesets %s",
        src,
        nrow(hm_data),
        ifelse(is_symmetric,
               " - symmetric similarity matrix\ncolor-coded by gene overlap (weighted by gene scores)",
               " - asymmetric similarity matrix\nscores reflect if columns are a supersets of rows")
        # at score=1, the superset/parent is shown in the column and child as row
        # High values imply respective columns contain most genes from geneset @ row
      )
    )
  }
}



#' Lollipop chart or barplot visualization of geneset enrichment testing results
#'
#' @examples \dontrun{
#'   # generate lollipop charts for each GO domain (CC/BP/MF), with geneset -log10
#'   # adjusted p-value on the x-axis and color-coding by geneset up/down-regulation
#'   plot_lollipop(
#'     result, output_dir = getwd(), plot_type = "lollipop", topn = 50,
#'     score_xaxis = "minlogp", score_color = "updown"
#'   )
#'
#'   # alternatively, as a barplot
#'   plot_lollipop(
#'     result, output_dir = getwd(), plot_type = "barplot", topn = 50,
#'     score_xaxis = "minlogp", score_color = "updown"
#'   )
#'
#'   # alternatively, color-code genesets by enrichment of significant genes using
#'   # parameter `score_color="oddsratio"` . See further `score_geneset_oddsratio`
#'   # function documentation for definition/computation of this score.
#'   plot_lollipop(
#'     result, output_dir = getwd(), plot_type = "lollipop", topn = 50,
#'     score_xaxis = "minlogp", score_color = "oddsratio"
#'   )
#' }
#' @param x results from function `test_genesets`
#' @param output_dir full path to output directory. Set to NA to directly show the figures instead of writing them to file
#' @param only_reduced only show the reduced/summarized set of significant genesets. This requires that you first applied the `reduce_genesets` function
#' @param plot_type Options: "barplot", "lollipop" (default)
#' @param show_pvalue boolean parameter that indicates whether adjusted p-values should be shown next to each data point
#' @param score_xaxis type of score to show on the x-axis. Options: "minlogp" to show -log10(adjusted pvalue), which is default. Use "oddsratio" to show the enrichment of significant genes.
#' For further details on this score and its computation, see the `score_geneset_oddsratio` function documentation.
#' Basically, the genesets in this plot are sorted by their proportion of foreground/significant genes
#' (and this ratio is standardized against the overall ratio of significant genes as to make this statistic comparable across analyses/datasets).
#' @param score_color analogous to `score_xaxis`, here you can specify the data used for color-coding the plot. Options: "minlogp", "oddsratio", "updown".
#' The former 2 options are the same as for `score_xaxis`, the latter enables color-coding by up-/down-regulation as encoded in the "score_type" column.
#' Note that this only works when geneset testing based on "effectsize" was performed and thus the "score_type" column has values "effectsize_up" or "effectsize_down" (encoding directionality).
#' Rows with other values are assumed NA and colored as grey.
#' @param score_color_limits defines the limits for the color scales. options; `score_color_limits = "source"` use the range of values per 'source' to compute colors scales (default). Set to "overall" in order to have a unified color scale across 'source'  (e.g. same color bar across GO_CC/GO_MF/GO_BP). Alternatively, provide a numeric vector of 2 values to manually define lower/upper-limits.
#' @param score_color_updown array of 2 strings that describe the colors for up- and down-regulation (used when `score_color` is set to "updown"). Default color-coding; up = red, down = blue. Use hex color codes, e.g. "#ff0000" for red.
#' @param max_ngenes only plot terms with less than N genes (quick way to get rid of large/unspecific terms)
#' @param topn topn terms to show after sorting genesets by p-value. For example, this makes it easy to plot the top10 GO terms. Set to NA to ignore this filter (default)
#' @param padj_cutoff adjusted pvalue cutoff for terms shown in the plot. If set to NA (default), all significant genesets are used (i.e. 'signif' column in the input geneset table)
#' @export
plot_lollipop = function(x, output_dir = NA, only_reduced = FALSE, plot_type = "lollipop", show_pvalue = FALSE, score_xaxis = "minlogp", score_color = ifelse(is.data.frame(x) && "score_type" %in% colnames(x) && is.character(x$score_type) && any(x$score_type %in% c("effectsize_up","effectsize_down")), "updown", "minlogp"), score_color_limits = "source", score_color_updown = c("#E57373", "#5E7ABC"), max_ngenes = NA, topn = NA, padj_cutoff = NA) {
  signif_and_reduced = ngenes = pvalue_adjust = pvalue = rank__ = name = score_oddsratio = score_clr = xvalue = clrvalue = NULL # fix invisible bindings R package NOTE
  check_dependency("ggplot2", "plot lollipop chart")
  x = validate_genesets(x)
  stopifnot("parameter output_dir must be a full path for the output files, this directory must already exists" =
              length(output_dir) == 1 && (is.na(output_dir) || (is.character(output_dir) && nchar(output_dir) > 4 && dir.exists(output_dir)) ) )
  stopifnot("parameter plot_type must be any of; 'barplot', 'lollipop'" = length(plot_type) == 1 && plot_type %in% c("barplot", "lollipop"))
  stopifnot("parameter score_xaxis must be any of; 'minlogp', 'oddsratio'" = length(score_xaxis) == 1 && score_xaxis %in% c("minlogp", "oddsratio"))
  stopifnot("parameter score_color must be any of; 'minlogp', 'oddsratio', 'updown'" = length(score_xaxis) == 1 && score_xaxis %in% c("minlogp", "oddsratio", "updown"))
  stopifnot("parameter score_color_limits a single string designating the color-scale scope (value must be either 'source' or 'overall'), or 2 numeric values that hardcode the color-scale limits" =
              (length(score_color_limits) == 1 && score_color_limits %in% c("source", "overall")) || (length(score_color_limits) == 2 && all(is.numeric(score_color_limits) & is.finite(score_color_limits)))  )
  stopifnot("parameter max_ngenes must be a single positive number (or NA to disable)" = length(max_ngenes) == 1 && (is.na(max_ngenes) || (is.numeric(max_ngenes) & is.finite(max_ngenes) & max_ngenes > 0)))
  stopifnot("parameter topn must be a single positive number (or NA to disable)" = length(topn) == 1 && (is.na(topn) || (is.numeric(topn) & is.finite(topn) & topn > 0)))
  stopifnot("parameter padj_cutoff must be a single positive number (or NA to disable)" = length(padj_cutoff) == 1 && (is.na(padj_cutoff) || (is.numeric(padj_cutoff) & is.finite(padj_cutoff) & padj_cutoff > 0)))
  stopifnot("parameter only_reduced must be a single boolean value" = length(only_reduced) == 1 && only_reduced %in% c(TRUE, FALSE))
  stopifnot("parameter only_reduced should not be used in conjunction with other filtering parameters like 'max_ngenes', 'topn', 'padj_cutoff'" = (only_reduced == FALSE) || (only_reduced && is.na(max_ngenes) && is.na(topn) && is.na(padj_cutoff)))

  if(score_xaxis == "oddsratio") {
    stopifnot("parameter score_xaxis was set to 'oddsratio' so we expect column 'score_oddsratio' in the provided data table (with finite numeric values)" = "score_oddsratio" %in% colnames(x) && all(is.finite(x$score_oddsratio)) )
  }
  if(score_color == "oddsratio") {
    stopifnot("parameter score_color was set to 'oddsratio' so we expect column 'score_oddsratio' in the provided data table (with finite numeric values)" = "score_oddsratio" %in% colnames(x) && all(is.finite(x$score_oddsratio)))
  }
  if(score_color == "updown") {
    stopifnot("parameter score_color was set to 'updown' so we expect column 'score_type' in the provided data table AND it should contain 'effectsize_up' and/or 'effectsize_down' values (i.e. when using GOAT with parameter score_type='effectsize')" = "score_type" %in% colnames(x) && all(is.character(x$score_type)) )
    if(!any(x$score_type %in% c("effectsize_up","effectsize_down"))) {
      cat("*** warning: parameter score_color was set to 'updown' but we did not find any 'effectsize_up' or 'effectsize_down' values in column 'score_type' and thus all elements will be color-coded in grey (NA). Try another color-coding, or use GOAT with parameter score_type='effectsize'\n")
    }
  }
  if(only_reduced) {
    stopifnot("'signif_and_reduced' column is not present in the data table, did you forget to use the reduce_genesets() function ?" = "signif_and_reduced" %in% colnames(x))
    x = x |> filter(signif_and_reduced %in% TRUE)
  }


  # optional user-provided filter criteria for the genesets to plot
  if(is.finite(max_ngenes)) {
    x = x |> filter(ngenes <= max_ngenes)
  }
  if(is.finite(padj_cutoff)) {
    x = x |> filter(pvalue_adjust <= padj_cutoff)
  } else {
    x = x |> filter(signif == TRUE)
  }
  if(is.finite(topn)) {
    x = x |>
      arrange(pvalue) |>
      group_by(source) |>
      mutate(rank__ = 1L:n()) |>
      ungroup() |>
      filter(rank__ <= topn) |>
      select(-rank__)
  }

  # return if there is no data to plot
  if(nrow(x) == 0) {
    cat("zero genesets match your filtering criteria, cannot generate plots\n")
    return()
  }

  # truncate geneset names
  label_nchar = 60
  x = x |> mutate(name = string_trunc_right(name, width = label_nchar))


  # define x-axis values
  xlab_title = color_title = "" # init
  x$xvalue = 0 # init column
  if(score_xaxis == "oddsratio") {
    xlab_title = "Significant gene odds ratio"
    x$xvalue = x$score_oddsratio
  } else {
    xlab_title = "-log10 adjusted p-value"
    x$xvalue = minlog10_fixzero(x$pvalue_adjust)
  }

  # analogous for color values
  clr_break_limits = NA # init val
  clr_break_limits_minmax = NA # sane default that defines minima/maxima
  x$clrvalue = 0 # init column
  if(score_color == "updown") {
    color_title = "Direction"
    x$clrvalue = "NA"
    x$clrvalue[x$score_type %in% "effectsize_up"] = "up"
    x$clrvalue[x$score_type %in% "effectsize_down"] = "down"
  } else if(score_color == "oddsratio") {
    color_title = "Odds ratio"
    x$clrvalue = x$score_oddsratio
    clr_break_limits_minmax = c(0.5, 2)
  } else {
    color_title = "-log10 p.adjust"
    clr_break_limits_minmax = c(0, 4)
    x$clrvalue = minlog10_fixzero(x$pvalue_adjust)
  }


  plotlist = list()

  # generate a plot for all 'source' elements
  for(src in sort(unique(x$source))) {
    # data subset for figure @ current source
    x_src = x |> filter(source == src)
    if(nrow(x_src) > 100) {
      cat(src, "has more than 100 significant genesets, not creating plot. You can run this function and explicitly request to plot only the topn genesets, see documentation for 'topn' parameter\n")
      next
    }

    # color scale
    if(length(score_color_limits) == 2) { # use hardcoded color scale limits if provided
      clr_break_limits = score_color_limits
    } else {
      if(score_color == "updown") {
        clr_break_limits = intersect(c("up", "down", "NA"), unique(x_src$clrvalue))
      } else {
        if(score_color_limits == "source") {
          clr_break_limits = c(min(clr_break_limits_minmax[1], min(x_src$clrvalue)), max(clr_break_limits_minmax[2], max(x_src$clrvalue)) )
        } else if(score_color_limits == "overall") {
          clr_break_limits = c(min(clr_break_limits_minmax[1], min(x$clrvalue)), max(clr_break_limits_minmax[2], max(x$clrvalue)) )
        } else {
          stop("unknown parameter for 'score_color_limits'")
        }
      }
    }

    # threshold colors to plot range
    if(length(clr_break_limits) == 2 && all(is.numeric(clr_break_limits))) {
      x_src$clrvalue[x_src$clrvalue < clr_break_limits[1]] = clr_break_limits[1]
      x_src$clrvalue[x_src$clrvalue > clr_break_limits[2]] = clr_break_limits[2]
    }

    # define y-axis plot order
    x_src = x_src |>
      arrange(desc(xvalue)) |>
      mutate(id = factor(id, levels = rev(unique(id))))

    # hack to always have same plot area for labels; pre-pad whitespace to fixed-width strings of N characters.
    # one could use stringi::stri_pad_left() or a wrapper thereof like stringr::str_pad(),
    # but instead we here use some naive vanilla R code to prevent an additional package dependency just for this.
    rows_pad = nchar(x_src$name) < label_nchar
    while(any(rows_pad)) {
      x_src$name[rows_pad] = paste0(" ", x_src$name[rows_pad])
      rows_pad = nchar(x_src$name) < label_nchar
    }

    # better implementation is available through labeling::extended() function, which is not included here to keep external dependencies to a minimum
    myxlim = c(0, max(pretty(c(0, max(1, max(x_src$xvalue))), n = 4)))
    size_breaks = c(10, 30, 100, 300, 1000)
    size_breaks = size_breaks[size_breaks <= max(x_src$ngenes)] # don't show "geneset size" legend elements that are larger than the actual data
    has_few_rows = nrow(x_src) < 10

    # create a very compact figure, reducing whitespace where possible
    p = NULL # declare val
    if(plot_type == "barplot") { # barplot
      p = ggplot2::ggplot(x_src, ggplot2::aes(x = xvalue, y = id, fill = clrvalue)) +
        ggplot2::geom_col(fill = "white") +
        ggplot2::geom_col()
    } else { # lollipop
      p = ggplot2::ggplot(x_src, ggplot2::aes(x = xvalue, y = id, size = ngenes, colour = clrvalue)) +
        ggplot2::geom_segment(ggplot2::aes(x=0, xend=xvalue, y=id, yend=id), colour = "darkgrey", size = 0.5) +
        ggplot2::geom_point(colour = "white", show.legend = FALSE) + # to plot alpha-colors as-is, first plot as white dot to overwrite panel + segment
        ggplot2::geom_point() +
        ggplot2::scale_size(range = c(2,6), limits = c(min(10, x_src$ngenes), max(1000, x_src$ngenes)), breaks = size_breaks, trans = "log10", name = "Count") +
        ggplot2::guides(size = ggplot2::guide_legend(keywidth = ggplot2::unit(0.1, "lines"), keyheight = ggplot2::unit(0.1, "lines"), order = 0))
    }

    p = p +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_x_continuous(limits = myxlim, expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_discrete(breaks = levels(x_src$id), labels = x_src$name[match(levels(x_src$id), x_src$id)]) +
      ggplot2::labs(x = xlab_title, y = "") +
      ggplot2::theme_bw()

    if(score_color == "updown") {
      p = p +
        # ifelse(plot_type == "barplot", color_title, "") # remove title for "up/down" to reduce required space
        ggplot2::scale_colour_manual(values = c("up"=score_color_updown[1], "down"=score_color_updown[2], "NA"="#7F7F7F"), breaks = clr_break_limits, name = color_title, aesthetics = c("colour", "fill")) +
        ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(size = 3)))
    } else {
      p = p + ggplot2::scale_colour_viridis_c(direction = 1, alpha = 0.7, limits = clr_break_limits, name = color_title, aesthetics = c("colour", "fill"))

      # deal with legends for very small plots
      # there seems to be a bug in ggplot2 that causes the legend keys to change based on the data / legend-element-count (first size then color legend, or vice versa)
      # so we hardcode the guides() here to always flow in the same ordering and use the 'order' parameter
      if(has_few_rows) {
        p = p + ggplot2::guides(
          colour = ggplot2::guide_colourbar(barwidth = ggplot2::unit(4, "lines"), order = 1, title.vjust = 1, label.theme = ggplot2::element_text(angle = -90, vjust = 0.5, size = 9))
        )
      } else {
        p = p + ggplot2::guides(colour = ggplot2::guide_colourbar(barheight = ggplot2::unit(4, "lines"), order = 1))
      }
    }

    p = p +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_line(size=0.5, colour = "black"),
        axis.ticks.y.left = ggplot2::element_blank(),
        axis.title.x.bottom = ggplot2::element_text(size = 11),
        legend.spacing.y = ggplot2::unit(0.1, "lines"), # reduce spacing between legend keys and items
        legend.title = ggplot2::element_text(size = 11),
        plot.margin = ggplot2::unit(c(0.5,1,0.5,0.5), "lines") # add some extra whitespace on the right side of the plot
      )

    if(has_few_rows) {
      p = p + ggplot2::theme(
        legend.position = "bottom",
        legend.margin = ggplot2::margin(2, 0, 0, 0), # reduce space between legend groups
      )
      if(plot_type == "lollipop") {
        p = p + ggplot2::theme(legend.justification = c(1, 0.5))
        print("test align")
      }
    }

    if(show_pvalue) {
      p = p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1e", pvalue_adjust)), size = 2.5, colour = "black", hjust = 0, nudge_x = max(x_src$xvalue)/15, show.legend = FALSE)
    }

    plotlist[[src]] = p

    # fiddle parameters to compute plot height;
    # <space needed for axis + vertical whitespace = independent of data>  +  <number of rows in the plot * X>  +  <optional padding if the legend is shown at the bottom>
    if(!is.na(output_dir)) {
      f = sprintf(
        "%s/goat_%s_%s%s%s%s%s_x=%s_clr=%s.pdf",
        output_dir,
        plot_type,
        gsub("[^a-zA-Z0-9 _-]+", "", src),
        ifelse(only_reduced, "_reduced", ""),
        ifelse(is.finite(max_ngenes), paste0("_", max_ngenes, "max_ngenes"), ""),
        ifelse(is.finite(topn), paste0("_", topn, "topn"), ""),
        ifelse(is.finite(padj_cutoff), paste0("_", padj_cutoff, "padj_cutoff"), ""),
        score_xaxis,
        score_color
      )

      if(file.exists(f)) {
        file.remove(f)
        if(file.exists(f)) {
          stop(paste("failed to overwrite existing file when creating a new plot:", f))
        }
      }

      ggplot2::ggsave(filename = f, plot = p, width = 7, height = 0.75 + nrow(x_src) * 0.14 + (has_few_rows * 0.7))
    }
  }

  if(is.na(output_dir)) {
    return(plotlist)
  }
}



#' plot geneset distance matrix as a network
#'
#' @param clusters result from `cluster_genesets`
#' @param src source property (e.g. "GO_CC")
#' @param show_clusters boolean value
#' @param show_text boolean value
#' @param topn_edges topN edges to retain per geneset (typically 5~8)
#' @param clr_default default color for the network, used only when `show_clusters` is set to `FALSE`
#' @export
plot_network = function(clusters, src, show_clusters = TRUE, show_text = FALSE, topn_edges = 5, clr_default = "#29b6f6") {
  source = cl = weight = genecount = cluster = label = NULL
  ## input validation
  # TODO

  ## prep data
  plot_title = src
  gs = clusters$genesets |> filter(source == src)
  mat = clusters$similarity[[src]]
  diag(mat) = 0

  ## from matrix to edgelist
  el = NULL
  for(i in 1:nrow(mat)) {
    i_cols = utils::head(order(mat[i,], decreasing = T), topn_edges)
    el = bind_rows(el, data.frame(from = rownames(mat)[i], to = colnames(mat)[i_cols], weight = mat[i,i_cols]) )
  }
  el = el |> filter(weight > 0.1)
  # transform weights to increase importance
  el$weight = el$weight^2

  ## create graph
  g = igraph::graph_from_data_frame(el)
  igraph::V(g)$genecount = clusters$genesets$ngenes[match(igraph::V(g)$name, clusters$genesets$id)]
  igraph::V(g)$label = clusters$genesets$name[match(igraph::V(g)$name, clusters$genesets$id)]

  ## layout
  set.seed(123)
  gl = igraph::layout_with_graphopt(g, charge = 0.05, spring.constant = 10)
  # gl = igraph::layout_with_fr(g, niter = 1000)

  ## cluster
  g_clust = igraph::cluster_optimal(g) # slow but reliable
  # g_clust = igraph::cluster_infomap(g)
  tmp = igraph::membership(g_clust)
  # cluster is assumed to be a continuous set of integers from 1 to N
  igraph::V(g)$cluster = as.integer(tmp[match(igraph::V(g)$name, names(tmp))])

  ## color
  df_clusters = data.frame(cl = unique(igraph::V(g)$cluster)) |> arrange(cl)
  df_clusters$clr = gg_color_hue(nrow(df_clusters))
  df_clusters$fill = lighten_color(df_clusters$clr, 0.25)

  p = ggraph::ggraph(g, layout = gl) +
    ggraph::geom_edge_link0(ggplot2::aes(width = I(0.5 + weight / 10),
                                         colour = I(grDevices::colorRampPalette(c("#ECEFF1", "#CFD8DC"))(10)[ cut(weight, breaks = (0:10)/10, include.lowest = T, labels=F) ]) ),
                            show.legend = F)

  if(show_clusters) {
    if(show_text) {
      p = p + ggraph::geom_node_point(ggplot2::aes(size = genecount, colour = I(df_clusters$fill[cluster])), show.legend = F)
    } else {
      p = p + ggraph::geom_node_point(ggplot2::aes(size = genecount, colour = I(df_clusters$clr[cluster]), fill = I(df_clusters$fill[cluster])), shape = 21, show.legend = F)
    }
  } else {
    p = p + ggraph::geom_node_point(ggplot2::aes(size = genecount), colour = clr_default, show.legend = F)
  }

  p = p +
    ggplot2::scale_size_continuous(trans = "log10", range = c(2, 6)) + # range = min/max node size
    ggplot2::guides(size = "none") +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_void() + # note; ggraph::theme_graph() is bugged
    ggplot2::theme(legend.position = "bottom",
                   legend.title = ggplot2::element_blank())

  if(show_text) {
    p = p + ggraph::geom_node_text(ggplot2::aes(label = label), size = 2)
  }

  return(p)
}
