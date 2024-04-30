
#' For each provided geneset, a volcano plot of all genelist log2fc and p-values with respective geneset constituents highlighted
#'
#' @examples \donttest{
#' # note; this example downloads data when first run, and typically takes ~60seconds
#'
#' # store the downloaded files in the following directory. Here, the temporary file
#' # directory is used. Alternatively, consider storing this data in a more permanent location.
#' # e.g. output_dir="~/data/goat" on unix systems or output_dir="C:/data/goat" on Windows
#' output_dir = tempdir()
#'
#' ## first run the default example from test_genesets() to obtain geneset results
#' datasets = download_goat_manuscript_data(output_dir)
#' genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`
#' genesets_asis = download_genesets_goatrepo(output_dir)
#' genesets_filtered = filter_genesets(genesets_asis, genelist)
#' result = test_genesets(genesets_filtered, genelist, method = "goat",
#'   score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)
#'
#' ## example 1; select top10 GO CC terms from the geneset testing results
#' result_subset = result |> filter(source == "GO_CC") |> arrange(pvalue) |> head(n = 10)
#' pdf(paste0(output_dir, "/volcano_CC_top10.pdf"), width = 4, height = 4)
#' plot_volcano(result_subset, genelist)
#' dev.off()
#'
#' ## example 2;, select small genesets that are significant and have
#' ## near-exclusive enrichment in either up up/down-regulated genes
#' # first, add geneset directionality scores to our results
#' result = score_geneset_directionality(result, genelist)
#' # next, subset the geneset result table
#' result_subset = result |>
#'   filter(signif & ngenes <= 50 & abs(score_directionality_rank) > 0.6) |>
#'   arrange(pvalue_adjust)
#' # finally, create plots. Note that the genelist contains a column 'symbol'
#' # which we use here to print labels for the topN genes per plotted geneset
#' pdf(paste0(output_dir, "/volcano_signif_ngenes50_directionality06.pdf"), width = 4, height = 4)
#' plot_volcano(result_subset, genelist, topn_labels = 10)
#' dev.off()
#' }
#' @param x a subset of the results from `test_genesets` function, see example
#' @param genelist input genelist, must contain columns 'gene', 'log2fc' and 'pvalue_adjust' (not! log transformed). If parameter topn_labels is provided, also include a character column 'symbol' that contains gene names/symbols/labels
#' @param plot if `TRUE`, will directly show the plots. if `FALSE`, returns a list of ggplot objects corresponding to rows in the input `result` parameter
#' @param topn_labels for how many genes that overlap between genelist and a geneset should we plot the gene symbol? This requires a column 'symbol' in the genelist parameter (default: 0)
#' @param color_default color for genes that are not part of a geneset (default: grey)
#' @param color_highlight color used to highlight geneset constituents (default: red)
#' @param color_label provided that topn_labels is set, this is the color of the text labels (default: black)
#' @param pointsize size of the dots, this parameter is passed to geom_point (default: 2)
#' @param pointalpha alpha of the dots, this parameter is passed to geom_point (default: 0.75)
#' @param labelsize provided that topn_labels is set, this is the text size (in pt) for the labels (default: 7)
#' @return if `plot==FALSE`, a list of ggplot2 objects. Otherwise, does not return any value
#' @export
plot_volcano = function(x, genelist, plot = TRUE, topn_labels = 0, color_default = "#B0B0B0", color_highlight = "#ef5350", color_label = "#000000", pointsize = 2, pointalpha = 0.75, labelsize = 7) {
  gene = log2fc = pvalue_adjust = highlight = show_label = symbol = NULL # fix invisible bindings R package NOTE
  check_dependency("ggplot2", "volcano plots")
  # input validation
  stopifnot("parameter x must be a result from the test_genesets function, a data.frame containing columns; source(character), id(character), name(character), genes(list), ngenes(integer), pvalue_adjust(numeric)" =
              is.data.frame(x) && nrow(x) > 0 && all(c("source", "id", "name", "genes", "ngenes", "pvalue_adjust") %in% colnames(x)) && is.list(x$genes) )
  stopifnot("parameter genelist must be the same genelist was was used as input for the test_genesets function AND contain columns required for creating volcano plots: gene (integer), log2fc(numeric) and pvalue_adjust(numeric)" =
              is.data.frame(genelist) && nrow(genelist) > 0 && all(c("gene", "log2fc", "pvalue_adjust") %in% colnames(genelist)) )
  stopifnot("parameter plot must be a single boolean value" = length(plot) == 1 && plot %in% c(TRUE, FALSE))
  stopifnot("parameter topn_labels must be a single non-negative integer value" = length(topn_labels) == 1 && is.numeric(topn_labels) && is.finite(topn_labels) && topn_labels >= 0)
  if(topn_labels != 0) {
    stopifnot("if parameter topn_labels is not zero, parameter topn_labels must contain a character column 'symbol'" = "symbol" %in% colnames(genelist) && is.character(genelist$symbol))
  }
  stopifnot("parameter color_default must be a single string representing a valid R color code" = length(color_default) == 1 && is.character(color_default) && nchar(color_default) > 2 && isvalid_color(color_default))
  stopifnot("parameter color_highlight must be a single string representing a valid R color code" = length(color_highlight) == 1 && is.character(color_highlight) && nchar(color_highlight) > 2 && isvalid_color(color_highlight))
  stopifnot("parameter color_label must be a single string representing a valid R color code" = length(color_label) == 1 && is.character(color_label) && nchar(color_label) > 2 && isvalid_color(color_label))
  stopifnot("parameter pointsize must be a single positive numeric value" = length(pointsize) == 1 && is.numeric(pointsize) && is.finite(pointsize) && pointsize > 0)
  stopifnot("parameter pointalpha must be a single non-negative numeric value" = length(pointalpha) == 1 && is.numeric(pointalpha) && is.finite(pointalpha) && pointalpha >= 0)
  stopifnot("parameter labelsize must be a single positive numeric value" = length(labelsize) == 1 && is.numeric(labelsize) && is.finite(labelsize) && labelsize > 0)
  topn_labels = as.integer(topn_labels)

  plotlist = list()
  for(i in 1:nrow(x)) {
    plotdata = genelist |>
      ungroup() |>
      select(gene, log2fc, pvalue_adjust, any_of("symbol")) |>
      mutate(
        show_label = FALSE,
        highlight = factor(gene %in% x$genes[[i]], levels = c("FALSE", "TRUE"))
      )
    # if labels are to be shown, sort the genelist by p-value and then flag the topn highlighted genes
    if(topn_labels != 0 && any(plotdata$highlight == TRUE)) {
      plotdata = plotdata |> arrange(pvalue_adjust)
      plotdata$show_label[utils::head(which(plotdata$highlight == TRUE), n = topn_labels)] = TRUE
    }
    # enforce plot ordering such that highlighted genes are on top
    plotdata = plotdata |> arrange(highlight, desc(pvalue_adjust))

    p = ggplot2::ggplot(plotdata, ggplot2::aes(log2fc, minlog10_fixzero(pvalue_adjust), colour = highlight)) +
      ggplot2::geom_point(size = pointsize, alpha = pointalpha, stroke = 0, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c("FALSE" = color_default, "TRUE" = color_highlight)) +
      ggplot2::labs(x = "log2 foldchange", y = "-log10 adjusted p-value",
                    title = sprintf("%s\n%s %s genes=%d adj.pvalue=%.2g", string_trunc_right(x$name[i], width = 70),
                                    x$source[i], ifelse(x$id[i] != x$name[i], x$id[i], ""), x$ngenes[i], x$pvalue_adjust[i])) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(size = 10),
        plot.title.position = "plot" # new param in ggplot 3.2.0 , can use because GOAT requires min version 3.3.0
      )

    if(topn_labels != 0 && any(plotdata$show_label == TRUE)) {
      p = p + ggplot2::geom_text(data = plotdata |> filter(show_label == TRUE), mapping = ggplot2::aes(label = symbol),
                                 size = labelsize / ggplot2::.pt, colour = color_label)
    }

    if(plot) {
      print(p)
    } else {
      plotlist[[i]] = p
    }
  }

  if(!plot) {
    return(plotlist)
  }
}
