
# 1) What do the gene-rank distributions of significant versus non-significant genesets look like in real-world data?
# - the rank values of genes from strongly enriched genesets look like a weak power-law (i.e. enrichment of top genes)
# - whereas genesets with no enrichment result in approximate uniform gene-rank distributions
#
# 2) ROC
# - generate 10k background genesets and 1k "genesets with some signal"  (we explore various patterns of enrichment)
# - compute geneset p-values with GOAT, GSEA and iDEA
# - perform ROC to test their sensitivity in separating the background genesets from those that should be enriched



#######################################################################################################################
############################################### compute & cache results ###############################################
#######################################################################################################################

require(pROC)
library(ggplot2)
library(goat) # also loads dplyr

# TODO: setup input/output paths
output_dir = getwd() # directory where output files should be stored
file_gene2go = "gene2go_2024-01-01.gz" # full path to NCBI gene2go file, e.g. previously downloaded from https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
file_goobo = "go_2024-01-01.obo" # full path go GO OBO file, e.g. previously downloaded from http://current.geneontology.org/ontology/go.obo
file_datasets = "goat_manuscript_datasets.rda" # full path to the datasets prepared earlier. This RData file is available in the "analyses" directory @ github
source("test_genesets_idea.R") # load the script that implements the iDEA function; update to the full path


# real-world genesets bundled with the GOAT R package that will be used for downstream generation of synthetic genesets
dataset_names_go_analysis = c(
  "Colameo 2021:RNA-seq:PMID34396684",
  "Sahadevan 2021:RNA-seq:PMID34021139",
  "Higginbotham 2020:mass-spec:PMID33087358",
  "Wingo 2020:mass-spec:PMID32424284"
)

# when the iDEA method is included, computations in this script will increase
# from approx. 1 hour to multiple weeks (!!) even on a fast workstation computer
include_idea = TRUE



goat_print_version()




############################################
############# helper functions #############
############################################

histogram_generank_distribution = function(gs, gl, gl_col_scoretype) {
  # 11 histogram breaks for 10 bins
  histbreaks = seq(from = 1, to = nrow(gl), length.out = 11)

  # for each geneset, compute histogram
  plotdata = NULL
  for(i in 1:nrow(gs)) {
    tmp = gl |> pull(!!gl_col_scoretype)
    i_hist_data = tmp[match(gs$genes[[i]], genelist$gene)]
    stopifnot(is.finite(i_hist_data))

    h = hist(i_hist_data, plot = F, breaks = histbreaks)

    plotdata = bind_rows(plotdata, tibble(
      id = gs$id[i],
      breaks = head(histbreaks, -1),
      density = h$counts / length(i_hist_data)
      # density = h$density
    ))
  }

  return(plotdata)
}



# sample k genes per genelist bin, hardcoded to enforce the gene rank distribution per generated geneset (over genelist bins)
generate_mock_genesets = function(list_genes_per_bin, genecount_per_bin, niter, geneset_id_template) {
  # prior to sampling, set.seed to ensure reproducible results
  set.seed(123)

  genesets = list()
  for(iter in seq_len(niter)) {
    # sample random genes
    iter_genes = NULL
    for(b in seq_along(list_genes_per_bin)) {
      iter_genes = c(iter_genes, sample(list_genes_per_bin[[b]], genecount_per_bin[b], replace = FALSE))
    }

    # create the mock geneset
    genesets[[length(genesets) + 1]] = tibble::tibble(
      source = "x",
      source_version = "x",
      id = paste0(geneset_id_template, iter),
      name = paste0(geneset_id_template, iter),
      genes = list(iter_genes),
      ngenes = length(iter_genes),
      ngenes_signif = 0L # not used in these analyses, enter a default value so the mock geneset is compatible with GOAT R package input geneset validations
    )
  }

  return(bind_rows(genesets))
}




############################################
######### check simulation designs #########
############################################

### visualize various distribution designs that we may use in above simulations
### i.e. these would go into `bins_foreground_all`

addpoints = function(x, col = 1, offset = 0) {
  stopifnot(sum(x) == 100) # double-check that we always use the same number of elements in this test
  lines(seq_along(x) + offset, x / sum(x), col = col)
  points(seq_along(x) + offset, x / sum(x), col = col, pch = 16, cex = 0.5)
}

plot(rep(0.1, 10), type = 'l', col = "grey", pch = 1, ylim = c(0.05, 0.2), xlab = 'gene rank "bin"', ylab = "Fraction")

addpoints(c(15, 12, 12, 11, 10, 9, 8, 8, 8, 7), col = "red", offset = -0.1) # strong enrichment of top-ranked genes, easy
addpoints(c(14, 12, 11, 10, 10, 9, 8, 8, 9, 9), col = "purple", offset = -0.2) # relatively easy to distinguish from uniform distribution because "many" genes are sampled from bins 1 and 2
addpoints(c(13, 11, 11, 10, 9, 10, 9, 9, 9, 9), col = "blue", offset = 0.1) # slightly more difficult, decreased number of genes from bins 1 and 2
addpoints(c(10, 13, 11, 10, 10, 9, 9, 9, 10, 9), col = "orange", offset = 0.2) # alternatively, enrichment pre top-genes
# most challenging; compared to uniform there will be a few more topN genes
# (i.e. will not yield typical "significant genesets", but methods should be able to distinguish this from the uniform distribution)
addpoints(c(12, 11, 11, 10, 10, 9, 9, 9, 9, 10), col = "green", offset = 0.3)




############################################
########### perform computations ###########
############################################

if(!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = F, recursive = T)
}

# load real-world OMICs datasets we prepared earlier from RData file
load(file_datasets)

# load GO genesets from gene2go file
genesets = load_genesets_go_fromfile(file_gene2go, file_goobo)

# iterate input genelists
dataset_results = list()
for(dataset_name in dataset_names_go_analysis) {
  # check if results are available from cache
  output_file_prefix = paste0(output_dir, "/", gsub("[^a-zA-Z0-9 _-]", " ", dataset_name))
  file_cache = paste0(output_file_prefix, "__roc_simulations_cache.rda")
  if(file.exists(file_cache)) {
    load(file_cache)
  } else { # no cache available, start computation
    cat("\nDATASET: ", dataset_name, "\n\n")
    time_start = Sys.time()
    gc(verbose = F, full = T)

    # genelist used in current iteration
    genelist = goat_manuscript_datasets[[dataset_name]]

    # add rank to genelist for downstream analyses
    tmp = goat_testgene_score_matrix(genelist, "effectsize")
    genelist$score_effectsize_up = tmp[,"effectsize_up"]
    genelist$score_effectsize_down = tmp[,"effectsize_down"]
    genelist = genelist |>
      arrange(desc(score_effectsize_up)) |> mutate(rank_effectsize_up = 1L:n()) |>
      arrange(desc(score_effectsize_down)) |> mutate(rank_effectsize_down = 1L:n())

    # filter genesets as per usual in GOAT workflow, removing genesets that have insufficient overlap with the input genelist
    genesets_filtered = filter_genesets(genesets, genelist, min_overlap = 10L, max_overlap = 1500L, max_overlap_fraction = 0.5, min_signif = NA, max_size = NA, dedupe = FALSE)

    # perform geneset testing using GSEA
    result = test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "effectsize")

    # now we can use these results to identify the topN and bottomN genesets and study their respective gene-rank distributions
    # here we are specifically using GSEA so the resulting geneset-to-generank-distribution profiles are not biased for GOAT


    ##################################################################################################
    ##################################################################################################


    ### iterate top/bottom-N genesets, compute histogram for each
    pct_top = 10
    pct_bottom = 10
    plotdata_empirical_distributions = list()

    for(myscoretype in c("effectsize_down", "effectsize_up")) {

      # genesets to consider; filter for direction-of-effectisze  +  geneset size
      tmp = result |> filter(score_type == myscoretype & ngenes >= 50 & ngenes <= 500)

      # foreground = top x% best genesets by respective p-value
      gs_foreground = tmp |> arrange(pvalue) |> slice_head(n = floor(nrow(tmp) * (pct_top/100)))

      # background = bottom x%
      gs_null = tmp |> arrange(desc(pvalue)) |> slice_head(n = floor(nrow(tmp) * (pct_bottom/100)))

      # compute histograms
      plotdata_foreground = histogram_generank_distribution(gs_foreground, genelist, paste0("rank_", myscoretype))
      plotdata_null = histogram_generank_distribution(gs_null, genelist, paste0("rank_", myscoretype))
      mtitle_foreground = sprintf("top %s%% genesets by '%s' and\n50 <= N <= 500 yields %d genesets", pct_top, myscoretype, nrow(gs_foreground))
      mtitle_null = sprintf("bottom %s%% genesets by '%s' and\n50 <= N <= 500 yields %d genesets", pct_bottom, myscoretype, nrow(gs_null))

      plotdata_empirical_distributions[[myscoretype]] = list(foreground = plotdata_foreground, null = plotdata_null,
                                                             foreground_title = mtitle_foreground, null_title = mtitle_null)
    }


    ##################################################################################################
    ##################################################################################################


    ### parameters
    niter_null = 10000L # number of mock genesets to generate for background set / null distribution
    niter_foreground = 1000L # number of mock genesets that should be distinct from the null
    bins_null = rep.int(10L, 10L) # generating null genesets, we draw 10 genes from 10 bins (index = bin, value = number of genes)
    # 5 templates for sampling genes to generate a foreground geneset
    bins_foreground_all = list(
      c(14, 12, 11, 10, 10, 9, 8, 8, 9, 9), # relatively easy to distinguish from uniform distribution because "many" genes are sampled from bins 1 and 2
      c(13, 11, 11, 10, 9, 10, 9, 9, 9, 9), # slightly more difficult, decreased number of genes from bins 1 and 2
      c(10, 13, 11, 10, 10, 9, 9, 9, 10, 9), # alternatively, enrichment in second bin / peak not in top-genes
      c(12, 11, 11, 10, 10, 9, 9, 9, 9, 10), # most challenging
      bins_null # null for suppl figure; ROC should be on diagonal
    )


    ### prepare plot data
    # split sorted input genelist into 10 bins (after sorting by effectsize)
    gl = genelist |>
      arrange(desc(score_effectsize_down)) |>
      mutate(
        rank = 1L:n(),
        bin = ggplot2::cut_interval(x = 1L:n(), n = length(bins_null), labels = F)
      )

    # cache gene identifiers per "bin" in a list for downstream computational efficiency
    tmp = gl |> select(gene, bin) |> tidyr::chop(cols = gene) |> arrange(bin)
    list_genes_per_bin = setNames(tmp$gene, tmp$bin)
    rm(tmp)


    ### generate mock datasets
    null_genesets = generate_mock_genesets(list_genes_per_bin, genecount_per_bin = bins_null, niter = niter_null, geneset_id_template = "gs_null_")

    ### for each template for mock genesets; generate synthetic genesets, then apply geneset testing with GOAT, GSEA, etc.
    roc_data = list()
    for(bin_index in seq_along(bins_foreground_all)) {
      bins_foreground = bins_foreground_all[[bin_index]]
      cat("\nDATASET: ", dataset_name, "testing geneset enrichment in bin", bin_index, "/", length(bins_foreground_all), "\n")

      # QC; ensure the geneset-frequency definitions align the number of bins and geneset-size between background and foreground
      stopifnot(length(bins_null) == length(bins_foreground) && sum(bins_null) == sum(bins_foreground))
      foreground_genesets = generate_mock_genesets(list_genes_per_bin, genecount_per_bin = bins_foreground, niter = niter_foreground, geneset_id_template = "gs_diff_")

      # compute geneset enrichment p-values for all mock genesets
      all_genesets = bind_rows(null_genesets, foreground_genesets)
      result_goat = test_genesets(all_genesets, gl, method = "goat", score_type = "effectsize")
      result_gsea = test_genesets(all_genesets, gl, method = "gsea", score_type = "effectsize")

      if(include_idea) {
        result_idea_default  = test_genesets(all_genesets, gl, method = "test_genesets_idea", idea_variant = "default", verbose = TRUE)
        result_idea_rescaled = test_genesets(all_genesets, gl, method = "test_genesets_idea", idea_variant = "rescale", verbose = TRUE)
        result_idea_beta_rescaled = test_genesets(all_genesets, gl, method = "test_genesets_idea", idea_variant = "beta_rescale", verbose = TRUE)
        # if("sd" %in% colnames(gl)) { # test/proof-of-concept; when log2fc SD are present in input data, try yet another iDEA variant
        #   result_idea_providedsd = test_genesets(all_genesets, gl, method = "test_genesets_idea", idea_variant = "genelist_sd", verbose = TRUE)
        #   roc_data[[length(roc_data) + 1]] = list(foreground_genesets = foreground_genesets, results = list(GOAT = result_goat, GSEA = result_gsea, iDEA = result_idea_default, "iDEA rescaled" = result_idea_rescaled, "iDEA beta rescaled" = result_idea_beta_rescaled, "iDEA provided SD" = result_idea_providedsd))
        # } else {
        roc_data[[length(roc_data) + 1]] = list(foreground_genesets = foreground_genesets, results = list(GOAT = result_goat, GSEA = result_gsea, iDEA = result_idea_default, "iDEA rescaled" = result_idea_rescaled, "iDEA beta rescaled" = result_idea_beta_rescaled))
        # }
      } else {
        roc_data[[length(roc_data) + 1]] = list(foreground_genesets = foreground_genesets, results = list(GOAT = result_goat, GSEA = result_gsea))
      }
    }

    ### finally, save to dataset-specific cache file
    save(dataset_name, output_file_prefix, gl, plotdata_empirical_distributions, bins_null, bins_foreground_all, null_genesets, roc_data,
         file = file_cache, compress = "xz")
    cat(sprintf("took; %.2f hours\n", as.numeric(difftime(Sys.time(), time_start, units = "hours"))))
  }



  ### append results for the current dataset to the final/aggregate results
  dataset_results[[dataset_name]] = list(
    dataset_name = dataset_name,
    output_file_prefix = output_file_prefix,
    genelist = gl,
    plotdata_empirical_distributions = plotdata_empirical_distributions,
    bins_null = bins_null,
    bins_foreground_all = bins_foreground_all,
    null_genesets = null_genesets,
    roc_data = roc_data
  )
} # end iteration over genelists


# finally, store all computation in a RData file
save(genesets, dataset_results, file = paste0(output_dir, "/roc_simulations.rda"), compress = "xz")













#######################################################################################################################
############################################### plot using cached data ################################################
#######################################################################################################################

require(pROC)
require(patchwork)
library(ggplot2)
library(goat) # also loads dplyr

# TODO: setup input/output paths
output_dir = getwd() # directory where output files should be stored
load(paste0(output_dir, "/roc_simulations.rda"))


############################################
############# helper functions #############
############################################

roc = function(l_result, legend_text_pt = 8, clr = c(GOAT = "#f50057", GSEA = "#304ffe", iDEA = "#C99D8E", "iDEA*" = "#B0BEC5", "iDEAb*" = "#D9A71E", "iDEAsd" = "#827717")) {
  ### compute ROCs, storing both the pROC object and the partial area under curve at 95% specificity
  # note that we order the genesets by their pvalue for ROC, i.e. multiple testing adjustment is ignored
  l_roc = list()
  l_roc_pauc = NULL
  for(n in names(l_result)) {
    l_result[[n]]$classification = c("foreground", "null")[grepl("null" ,l_result[[n]]$id) + 1]
    tmp = pROC::roc(l_result[[n]]$classification, abs(l_result[[n]]$pvalue), levels = c("null", "foreground"), direction=">", partial.auc=c(1, 0.95), partial.auc.correct=F, partial.auc.focus="specificity")
    l_roc_pauc = c(l_roc_pauc, as.numeric(tmp$auc))
    l_roc[[n]] = pROC::roc(l_result[[n]]$classification, abs(l_result[[n]]$pvalue), levels = c("null", "foreground"), direction=">")
  }
  l_roc_legend = sprintf("%s %.2f%%", names(l_roc), l_roc_pauc * 100)


  ### ggplot
  p = pROC::ggroc(l_roc, aes = "colour") +
    geom_vline(xintercept = 0.95, color = "darkgrey") +
    geom_segment(mapping = aes(x = 1, y = 0, xend = 0, yend = 1), data = data.frame(), color = "black", linetype = "dashed", size = 1) +
    scale_colour_manual(values = clr,
                        breaks = names(l_result),
                        labels = l_roc_legend) +
    coord_cartesian(xlim = 1:0, ylim = 0:1) +
    # scale_x_continuous(limits = 0:1, expand = expansion(mult = c(0.05, 0.05))) +
    # scale_y_continuous(limits = 0:1, expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Specificity", y = "Sensitivity") +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.position = "none"
    )

  # custom legend by manually coding ggplot annotations, allowing use to use colored text (reduce legend size because we no longer need the lines)
  # 8pt text, white background, no border
  for(i in seq_along(l_result)) {
    i_clr = clr[names(l_result)[i]]
    if(is.na(i_clr)) i_clr = "grey"
    p = p + annotate(geom = "label", x = 0.55, y = 0.01 + 0.1 * (i-1), label = l_roc_legend[i], size = legend_text_pt/.pt, hjust = 0, color = i_clr, fill = NA, label.size = NA) # , fill = "white"
  }

  return(p)
}


custom_density_plot = function(g, add_normal = TRUE) {
  p = g +
    geom_vline(xintercept = 0, colour = "darkgrey") +
    geom_line(stat = "density", adjust = 0.5)
  if(add_normal) {
    p = p +
      geom_line(data = data.frame(x = seq(from = -5, to = 5, by = 0.1)), aes(x, dnorm(x)), inherit.aes = FALSE, color = "black")
  }
  p +
    guides(colour = guide_legend(nrow = 2)) +
    labs(y = "Density") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      plot.background = element_blank(),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11)
    )
}




############################################
########### generate ROC figures ###########
############################################

for(dataset_name in names(dataset_results)) {
  x = dataset_results[[dataset_name]]

  plotlist_distr = plotlist_roc = plotlist_pval_dens = list()
  plotlist_pval_dens__npanel = NULL
  for(index_bins in seq_along(x$bins_foreground_all)) {
    bins_foreground = x$bins_foreground_all[[index_bins]]

    ### visualize the mock geneset distributions (i.e. how many genes were sampled from each bin in the genelist)

    gl_bins = x$genelist |> distinct(bin, .keep_all = T) |> select(bin, rank)
    tib_plot = tibble(bin = c(seq_along(x$bins_null), seq_along(bins_foreground)),
                      count = c(x$bins_null, bins_foreground),
                      frac = c(x$bins_null/sum(x$bins_null), bins_foreground/sum(bins_foreground)),
                      type = c(rep("null", length(x$bins_null)), rep("foreground", length(x$bins_null))) ) |>
      mutate(generank = gl_bins$rank[match(bin, gl_bins$bin)])

    p_distribution = ggplot(tib_plot, aes(x = generank, y = frac, colour = type, shape = type)) +
      geom_line(show.legend = F) +
      geom_point(size = 2, show.legend = F) +
      scale_colour_manual(values = c(null = "#666666", foreground = "#4caf50")) +
      scale_shape_manual(values = c(null = 16, foreground = 17)) +
      scale_x_continuous(n.breaks = 4) +
      coord_cartesian(ylim = c(0.05, 0.18)) +
      labs(x = "Gene rank", y = "Fraction") +
      theme_bw() +
      theme(
        plot.background = element_blank(),
        axis.text = element_text(size = 9),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9)
      )

    ### ROC based on plain geneset p-values (sans multiple testing adjustment)
    l = x$roc_data[[index_bins]]
    plotdata = list(
      GOAT = l$results$GOAT |> select(id, pvalue),
      GSEA = l$results$GSEA |> select(id, pvalue)
    )
    # for iDEA, only use louis pvalues in ROC
    if("iDEA" %in% names(l$results)) {
      plotdata[["iDEA"]] = l$results$iDEA |> select(id, pvalue = pvalue_louis_idea)
    }
    if("iDEA rescaled" %in% names(l$results)) {
      plotdata[["iDEA*"]] = l$results$`iDEA rescaled` |> select(id, pvalue = pvalue_louis_idea)
    }
    if("iDEA beta rescaled" %in% names(l$results)) {
      plotdata[["iDEAb*"]] = l$results$`iDEA beta rescaled` |> select(id, pvalue = pvalue_louis_idea)
    }
    if("iDEA provided SD" %in% names(l$results)) {
      plotdata[["iDEAsd"]] = l$results$`iDEA provided SD` |> select(id, pvalue = pvalue_louis_idea)
    }

    p_roc = roc(plotdata)



    ### p-value distribution plots
    l = x$roc_data[[index_bins]]
    plotdata = bind_rows(
      l$results$GOAT |> select(id, pvalue) |> mutate(method = "GOAT"),
      l$results$GSEA |> select(id, pvalue) |> mutate(method = "GSEA")
    )
    if("iDEA" %in% names(l$results)) {
      plotdata = bind_rows(
        plotdata,
        l$results$iDEA |> select(id, pvalue = pvalue_idea) |> mutate(method = "iDEA as-is"),
        l$results$iDEA |> select(id, pvalue = pvalue_louis_idea) |> mutate(method = "iDEA louis")
      )
    }
    if("iDEA rescaled" %in% names(l$results)) {
      plotdata = bind_rows(
        plotdata,
        l$results$`iDEA rescaled` |> select(id, pvalue = pvalue_idea) |> mutate(method = "iDEA* as-is"),
        l$results$`iDEA rescaled` |> select(id, pvalue = pvalue_louis_idea) |> mutate(method = "iDEA* louis")
      )
    }
    if("iDEA beta rescaled" %in% names(l$results)) {
      plotdata = bind_rows(
        plotdata,
        l$results$`iDEA beta rescaled` |> select(id, pvalue = pvalue_idea) |> mutate(method = "iDEAb* as-is"),
        l$results$`iDEA beta rescaled` |> select(id, pvalue = pvalue_louis_idea) |> mutate(method = "iDEAb* louis")
      )
    }
    if("iDEA provided SD" %in% names(l$results)) {
      plotdata = bind_rows(
        plotdata,
        l$results$`iDEA provided SD` |> select(id, pvalue = pvalue_idea) |> mutate(method = "iDEAsd as-is"),
        l$results$`iDEA provided SD` |> select(id, pvalue = pvalue_louis_idea) |> mutate(method = "iDEAsd louis")
      )
    }

    plotdata = plotdata |> mutate(method = factor(method, levels = unique(method)))
    p_pvaldens = ggplot(plotdata, aes(x = pvalue, colour = grepl("null", id))) +
      geom_density(show.legend = F) +
      scale_colour_manual(values = c("TRUE" = "#666666", "FALSE" = "#4caf50")) + # TRUE = null, FALSE = foreground
      facet_wrap(.~method, scales = "free_y", ncol = 1) +
      labs(x = "p-value", y = "Density") +
      theme_bw() +
      theme(
        plot.background = element_blank(),
        axis.text = element_text(size = 9),
        strip.text = element_text(size = 11),
        axis.title = element_text(size = 11)
      )


    ### finally, store data
    plotlist_distr[[length(plotlist_distr) + 1]] = p_distribution
    plotlist_roc[[length(plotlist_roc) + 1]] = p_distribution
    plotlist_roc[[length(plotlist_roc) + 1]] = p_roc
    plotlist_pval_dens[[length(plotlist_pval_dens) + 1]] = p_pvaldens
    plotlist_pval_dens__npanel = c(plotlist_pval_dens__npanel, length(levels(plotdata$method)))
  }


  ### empirical data

  plot_empirical = function(x) {
    ggplot(x, aes(x = breaks, y = density)) +
      geom_line(aes(group = id), alpha = 0.25, colour = "#777777") +
      geom_smooth(formula = y ~ x, method = "loess", se = TRUE, span = 0.5, level = 0.95, colour = "#ff1744", fill = "#40c4ff", size = 1.25) +
      scale_x_continuous(n.breaks = 4) +
      labs(x = "Gene rank", y = "Fraction") +
      theme_bw() +
      theme(
        plot.background = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 11)
      )
  }

  p_empirical_foreground = plot_empirical(x$plotdata_empirical_distributions$effectsize_down$foreground)
  p_empirical_null = plot_empirical(x$plotdata_empirical_distributions$effectsize_down$null)


  ### construct combined plot

  # only the first 2 templates, i.e. the easy ones
  plotlist_main = list()
  plotlist_main[[1]] = p_empirical_foreground
  plotlist_main[[2]] = p_empirical_null
  for(i in c(1:4, length(plotlist_roc) - 1, length(plotlist_roc))) plotlist_main[[length(plotlist_main) + 1]] = plotlist_roc[[i]]

  # suppl. data; show all templates
  plotlist_suppl = list()
  plotlist_suppl[[1]] = p_empirical_foreground
  plotlist_suppl[[2]] = p_empirical_null
  for(i in seq_along(plotlist_roc)) plotlist_suppl[[length(plotlist_suppl) + 1]] = plotlist_roc[[i]]


  p_combined_main = patchwork::wrap_plots(plotlist_main, nrow = 2, byrow = FALSE) +
    patchwork::plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 13))

  p_combined_suppl = patchwork::wrap_plots(plotlist_suppl, nrow = 4, byrow = FALSE) +
    patchwork::plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 13))


  ggsave(paste0(x$output_file_prefix, "__roc_simulations.pdf"), plot = p_combined_main, width = 8.5, height = 4)
  ggsave(paste0(x$output_file_prefix, "__roc_simulations__suppl.pdf"), plot = p_combined_suppl, width = 6.5, height = 7.5)



  ### p-value distribution plots

  plotlist = list()
  N = length(plotlist_distr)
  # N = length(plotlist_distr) - 1 # optionally, skip the null distribution plot to save some horizontal space
  for(i in 1:N) {
    if(i == 1) plotlist[[length(plotlist) + 1]] = plotlist_distr[[i]] + theme(plot.margin = grid::unit(c(5.5, 1, 5.5, 5.5), "pt"), axis.text = element_text(size = 7))
    else plotlist[[length(plotlist) + 1]] = plotlist_distr[[i]] + labs(y = "") + theme(plot.margin = grid::unit(c(5.5, 1, 5.5, 0.5), "pt"), axis.text = element_text(size = 7))
  }
  for(i in 1:N) {
    if(i == 1) plotlist[[length(plotlist) + 1]] = plotlist_pval_dens[[i]] + theme(plot.margin = grid::unit(c(5.5, 1, 5.5, 5.5), "pt"), axis.text = element_text(size = 7))
    else plotlist[[length(plotlist) + 1]] = plotlist_pval_dens[[i]] + labs(y = "") + theme(plot.margin = grid::unit(c(5.5, 1, 5.5, 0.5), "pt"), axis.text = element_text(size = 7))
  }

  p_combined_pvaldens = patchwork::wrap_plots(plotlist, nrow = 2, ncol = N, byrow = TRUE, heights = c(1, max(plotlist_pval_dens__npanel))) +
    patchwork::plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 13, hjust = 1))

  ggsave(paste0(x$output_file_prefix, "__roc_simulations__pvalue-density.pdf"), plot = p_combined_pvaldens, width = 8.5,
         height = ifelse(max(plotlist_pval_dens__npanel) == 2, 4.5, ifelse(max(plotlist_pval_dens__npanel) <= 4, 6, ifelse(max(plotlist_pval_dens__npanel) <= 6, 7.5, 9))) )

} # end iteration over genelists







############################################
#### compare iDEA input between datasets ###
############################################

idea_input_data = NULL
for(dataset_name in names(dataset_results)) {
  genelist = goat_manuscript_datasets[[dataset_name]]

  ### below code snippet is copy/pasted from our iDEA implementation @ test_genesets_idea.R
  # here we basically follow the iDEA tutorial @ https://xzhoulab.github.io/iDEA/, but with minor adaptions for edge-cases that cause zero or infinite se_beta estimates
  genelist_idea = data.frame(beta = genelist$log2fc, beta_var = NA, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
  genelist_idea$zscore = stats::qnorm(genelist_idea$pvalue / 2.0, lower.tail=FALSE)
  genelist_idea$se_beta = abs(genelist_idea$beta / genelist_idea$zscore) # approximate the standard error of beta (log2fc)
  genelist_idea$se_beta[!is.finite(genelist_idea$se_beta)] = 10 # update/fix; default high value where beta_se is missing, yields beta_var = 100
  genelist_idea$se_beta[genelist_idea$se_beta < 10^-6] = 10^-6 # guard against "zero standard error" (e.g. at log2fc=0; 0/x=0) to ensure downstream iDEA code doesn't suffer from division-by-zero issues
  genelist_idea$beta_var = genelist_idea$se_beta^2

  idea_input_data = bind_rows(idea_input_data, genelist_idea |> mutate(dataset = dataset_name))
}

plotdata = idea_input_data |> mutate(dataset_label = sub(":", ", ", sub(":PMID.*", "", dataset)))
# input data, issue 1; beta/log2fc is not approximately normally distributed. Higginbotham-MS = outlier here, has second-most NA values @ simulations
p1 = custom_density_plot(ggplot(plotdata, aes(beta, colour = dataset_label)), add_normal = F) +
  coord_cartesian(c(-2.5, 2.5))
# input data, issue 2: beta/beta_var differ greatly between datasets. Wingo is a clear outlier
p2 = custom_density_plot(ggplot(plotdata, aes(beta / beta_var, colour = dataset_label)), add_normal = F) +
  coord_cartesian(c(-50, 50))

p_idea_betavar_distr = patchwork::wrap_plots(p1, p2, nrow = 1, ncol = 2, byrow = TRUE, guides = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 13), legend.position = "bottom")

ggsave(paste0(output_dir, "/idea_input_data.pdf"), plot = p_idea_betavar_distr, width = 5, height = 3.5)


