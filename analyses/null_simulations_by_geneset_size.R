
library(ggplot2)
library(goat) # also loads dplyr

# TODO: setup input/output paths
output_dir = getwd() # directory where output files should be stored
file_gene2go = "gene2go_2024-01-01.gz" # full path to NCBI gene2go file, e.g. previously downloaded from https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
file_goobo = "go_2024-01-01.obo" # full path go GO OBO file, e.g. previously downloaded from http://current.geneontology.org/ontology/go.obo
file_datasets = "goat_manuscript_datasets.rda" # full path to the datasets prepared earlier. This RData file is available in the "analyses" directory @ github
source("test_genesets_idea.R") # load the script that implements the iDEA function; update to the full path


# hardcoded color-codings
clr = tibble::tibble(
  method = c("hypergeometric", "fisherexact_ease", "idea_default_louis", "idea_rescale_louis", "idea_beta_rescale_louis",
             "gsea_pvalue", "gsea_effectsize", "gsea_pvalue_gseaParam2", "gsea_effectsize_gseaParam2",
             "goat_pvalue", "goat_effectsize", "goat_pvalue_effectsize", "goat_effectsize_up", "goat_effectsize_down", "goat_effectsize_abs"),
  label = c("hypergeometric test", "Fisher-exact EASE", "iDEA", "iDEA*", "iDEAb*",
            "GSEA p-value", "GSEA effect size", "GSEA p-value param=2", "GSEA effectsize param=2",
            "GOAT p-value", "GOAT effect size", "GOAT p-value & effect size", "GOAT effect size up", "GOAT effect size down", "GOAT effect size abs"),
  clr = c("#69f0ae", "#ffd600", "#C99D8E", "#B0BEC5", "#D9A71E",
          "#42a5f5", "#304ffe", "#5c6bc0", "#5e35b1",
          "#ff7043", "#f50057", "#d50000", "#f50057", "#f50057", "#f50057")
)


goat_print_version()




##########################################################
### GENERATE SYNTHETIC NULL GENESETS + APPLY GOAT/GSEA ###
##########################################################


benchmark_bin_results = function(x) {
  y = x |>
    # 'best' pvalue on top
    arrange(pvalue) |>
    # group by geneset size ('ngenes' column)
    group_by(ngenes) |>
    mutate(
      # frac = fraction of genesets seen at row N  (relative within grouping) = expected p
      frac = (1:n()) / n(),
      pvalue_minlog10 = minlog10_fixzero(pvalue),
      frac_minlog10 = minlog10_fixzero(frac),
      # bin the data (within grouping) within 0.05 distance on -log10 scale  (so we don't have to plot 200k datapoints later on)
      bin = cut(frac_minlog10, breaks = seq(from=0, to=ceiling(max(frac_minlog10)), by = 0.05), include.lowest = T, labels = F)
    )

  # now that we have observed (pvalue) and expected (frac) -log10 values, we can quantify the error
  # here we compute the RMSE within "ngenes" grouping
  # (because function input are the geneset testing results from 1 algorithm, and here the table is grouped by "ngenes", this computes the respective error per geneset size)
  err_ngenes = y |> summarise(rmse_ngenes = sqrt(mean((frac_minlog10 - pvalue_minlog10)^2)), .groups = "drop")
  # analogous for overall data (i.e. over all genesets)
  err_overall = sqrt(mean((y$frac_minlog10 - y$pvalue_minlog10)^2))
  # err_overall = y |> ungroup() |> summarise(rmse_overall = sqrt(mean((frac_minlog10 - pvalue_minlog10)^2)), .groups = "drop")

  # summary statistics, reducing the huge table to only information we will plot (1 value per ngenes*bin)
  y = y |>
    group_by(ngenes, bin) |>
    summarise(
      count_genesets = n(),
      frac_minlog10 = mean(frac_minlog10),
      pvalue_minlog10 = mean(pvalue_minlog10),
      .groups = "drop"
    )

  y |>
    left_join(err_ngenes, by = "ngenes") |>
    mutate(rmse_overall = err_overall)
}



file_cache = paste0(output_dir, "/null_simulations.rda")
if(file.exists(file_cache)) {
  load(file_cache)
} else {
  # load real-world OMICs datasets we prepared earlier from RData file
  load(file_datasets)

  # setup input data; GO genesets  &  real-world dataset that is to be tested against randomized genesets
  genesets = load_genesets_go_fromfile(file_gene2go, file_goobo)

  # setup input data; we'll sample from a real-world genelist
  genelist = goat_manuscript_datasets$`Sahadevan 2021:RNA-seq:PMID34021139`

  # identifiers for the random genesets we'll generate. typically some vector of length 100k ~ 500k
  # the length of this vector defines the number of random genesets we'll generate
  idseq = 1L:200000L # 200k geneset permutations

  # ensure our results are reproducible
  set.seed(123)

  tib_plot_all = NULL
  for(genelist_size in c(250L, 1000L, 4000L, 12000L)) {
    stopifnot(genelist_size < nrow(genelist))

    cat("genelist_size =", genelist_size, "\n")

    ### generate the mock data

    # subset N random genes from the genelist
    mock_genelist = genelist |> slice(sample(seq_len(nrow(genelist)), size = genelist_size))

    # generate many random genesets of size k
    genesets_mock = NULL
    for(mock_size in c(10L, 20L, 50L, 100L, 200L, 1000L)) {
      # limit geneset sizes to at most half the genelist length (larger = go next)
      if(mock_size * 2 > genelist_size) {
        next
      }

      # generate mock genesets; draw k random gene IDs from the genelist. Importantly, sample without replacement
      genesets_mock = bind_rows(
        genesets_mock,
        tmp = tibble::tibble(
          source = "x",
          source_version = "x",
          id = paste0("gs_", mock_size, "_", idseq),
          name = id,
          genes = sapply(idseq, function(i) sample(mock_genelist$gene, size = mock_size, replace = FALSE), simplify = FALSE, USE.NAMES = FALSE),
          ngenes = mock_size,
          ngenes_signif = 0L
        )
      )
    }


    ### APPLY TEST_GENESETS
    # separate statements, instead of 1 call to bind_rows, so we don't have to hold intermediate results across all algorithms in memory (i.e. R can clean up unused RAM in between)

    # GOAT
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "pvalue")) |>
                               mutate(genelist_size = genelist_size, method_name = "GOAT p-value"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize")) |>
                               mutate(genelist_size = genelist_size, method_name = "GOAT effect size"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_up")) |>
                               mutate(genelist_size = genelist_size, method_name = "GOAT effect size up"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_down")) |>
                               mutate(genelist_size = genelist_size, method_name = "GOAT effect size down"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_abs")) |>
                               mutate(genelist_size = genelist_size, method_name = "GOAT effect size abs"))
    # GSEA
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 1000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=1000"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 1000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA effect size\nniter=1000"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 10000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=10000"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 10000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA effect size\nniter=10000"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 50000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=50000"))
    tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 50000L)) |>
                               mutate(genelist_size = genelist_size, method_name = "GSEA effect size\nniter=50000"))

    ### debug; apply GOAT with specific parameters, then create QQplot that splits results per geneneset size/bin
    ### debug; color-code by the number of observed datapoints in each size/bin -->> be mindful that this highlights/emphasizes one-hit-wonders / extreme values
    # x = test_genesets_goat_precomputed(genesets = genesets_mock, genelist = mock_genelist, score_type = "pvalue")
    # y = benchmark_bin_results(x)
    # p = ggplot(y, aes(frac_minlog10, pvalue_minlog10, colour = pmin(count_genesets, 500))) +
    #   geom_point(size = 0.3) +
    #   geom_abline(intercept = 0, slope = 1, colour = "red", size = 0.5) +
    #   facet_wrap(.~ngenes) +
    #   scale_colour_continuous(name = "# data points", trans = "log2", high = "#000000", low = "#BBBBBB") +
    #   labs(title = paste("genelist N=", genelist_size)) +
    #   theme_bw() +
    #   theme(legend.position = "bottom")
    # print(p)
  }


  ### finally, store computations in RData file
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  save(tib_plot_all, file = file_cache, compress = TRUE)
}






##########################################################
##################### CREATE FIGURES #####################
##########################################################


# Plot only geneset bins where we have at least 25 datapoints to prevent extreme values,
# that are observed only once or twice out of 200 thousand (!) random draws,
# from showing up in the plot and causing strong bias in interpretation
# Outliers not supported by N observations are not indicative of an algorithm's performance.
# (note that with different set.seed, those 1-out-of-200k extremes will vary while the rest should not due to averaging across bootstrap iterations)
myplot = function(x, min_datapoints = 25, show_rmse = TRUE) {
  y = x |>
    filter(count_genesets >= min_datapoints) |>
    # update method names
    mutate(
      method_label = gsub("_+", " ", method_name),
      method_label = factor(method_label, levels = unique(method_label))
    ) |>
    arrange(desc(ngenes)) |>
    mutate(ngenes = factor(ngenes, unique(ngenes))) |>
    arrange(genelist_size) |>
    mutate(
      genelist_size_label = paste0("N = ", genelist_size),
      genelist_size_label = factor(genelist_size_label, unique(genelist_size_label))
    )

  p = ggplot(y, aes(frac_minlog10, pvalue_minlog10, colour = ngenes)) +
    geom_line(size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", size = 0.4) +
    facet_wrap(method_label ~ genelist_size_label, ncol = n_distinct(y$method_label), dir = "v") +
    guides(colour = guide_legend(title = "gene set size", nrow = 1, reverse = T, override.aes = list(size = 1))) +
    scale_colour_discrete(direction = -1) +
    labs(x = "Expected -log10 p-value", y = "Observed -log10 p-value") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 8, margin = margin(2,1,2,1)),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )

  if(show_rmse) {
    # from input data, grab the error for each unique genelist*method pair
    err = y |>
      distinct(genelist_size, method_name, .keep_all = TRUE) |>
      mutate(err_label = sprintf("RMSE: %.4f", rmse_overall))
    # add annotation to the plot on the top-left
    p = p +
      geom_label(mapping = aes(x = -0.05, y = 3.3, label = err_label, colour = I("black")), data = err, size = 7/ggplot2::.pt, hjust = 0, vjust = 1, fill = "white", label.size = NA, show.legend = FALSE)
  }

  return(p)
}



# GOAT pvalue, GOAT effectsize, GSEA 1kiter pvalue, GSEA 1kiter effectsize, GSEA 50kiter pvalue, GSEA 50kiter effectsize
plot_names = c("GOAT p-value", "GOAT effect size", "GSEA p-value\nniter=1000", "GSEA effect size\nniter=1000", "GSEA p-value\nniter=50000", "GSEA effect size\nniter=50000")
p = myplot(tib_plot_all |> filter(method_name %in% plot_names) |> arrange(match(method_name, plot_names)) )
p = p + facet_grid(vars(genelist_size_label), vars(method_label))
print(p)
ggsave(p, file = paste0(output_dir, "/null_simulations_bins__figure2.pdf"), width = 6.5, height = 5)

# all GOAT score types
p = myplot(tib_plot_all |> filter(grepl("GOAT", method_name)) )
print(p)
ggsave(p, file = paste0(output_dir, "/null_simulations_bins__GOAT-all-scoretypes.pdf"), width = 7, height = 7)

# GSEA 1k/10k/50k iter compared
p = myplot(tib_plot_all |> filter(grepl("GSEA", method_name)) )
print(p)
ggsave(p, file = paste0(output_dir, "/null_simulations_bins__GSEA-compare-niter.pdf"), width = 7, height = 7)

# p = myplot(tib_plot_all |> filter(method_name == "GOAT p-value")) +
#   facet_wrap(method_label ~ genelist_size_label, nrow = 1)
# svglite::svglite(filename = "C:/temp/null_simulations.svg", width = 8, height = 3)
# print(p)

# only GOAT p-value and effectsize
p = myplot(tib_plot_all |> filter(method_name == "GOAT p-value" | method_name == "GOAT effectsize")) +
  facet_wrap(method_label ~ genelist_size_label, nrow = 2)
print(p)
ggsave(p, file = paste0(output_dir, "/null_simulations_bins__GOAT-pvalue-effectsize.pdf"), width = 5, height = 3.5)









######### RMSE values for all method*genelist*ngenes combinations (same methods as null simulations main figure)

plot_names = c("GOAT p-value", "GOAT effect size", "GSEA p-value\nniter=1000", "GSEA effect size\nniter=1000", "GSEA p-value\nniter=50000", "GSEA effect size\nniter=50000")
x = tib_plot_all |> filter(method_name %in% plot_names) |> arrange(match(method_name, plot_names))

y = x |>
  distinct(method_name, genelist_size, ngenes, .keep_all = T) |>
  select(method_name, genelist_size, ngenes, rmse_ngenes) |>
  # below labels and factors are analogous to main figure plot/code
  arrange(desc(ngenes)) |>
  mutate(ngenes = factor(ngenes, unique(ngenes))) |>
  arrange(genelist_size) |>
  mutate(
    genelist_size_label = paste0("N = ", genelist_size),
    genelist_size_label = factor(genelist_size_label, unique(genelist_size_label)),
    classification = factor(ifelse(grepl("value", method_name), "p-value", "effect size"), levels = c("p-value", "effect size"))
  )

data_legend = y |>
  distinct(method_name) |>
  mutate(
    method_clr_label = gsub("\n.*", "", method_name),
    method_clr = clr$clr[match(gsub("_", " ", method_clr_label), clr$label)],
    method_shape = ifelse(grepl("niter=1000", method_name), 15, ifelse(grepl("niter=50000", method_name), 17, 16)),
    method_label = factor(method_name, levels = unique(method_name)) # follow order defined upstream
  ) |>
  arrange(method_label)
stopifnot(!is.na(data_legend$method_clr)) # double-check that the method name lookup succeeded

p = ggplot(y, aes(rmse_ngenes, ngenes, colour = method_name, shape = method_name)) +
  geom_point(size = 2, alpha = 0.66) +
  scale_colour_manual(name = "method", values = tibble::deframe(data_legend |> select(method_name, method_clr)), breaks = data_legend$method_name, labels = data_legend$method_name) +
  scale_shape_manual(name = "method", values = tibble::deframe(data_legend |> select(method_name, method_shape)), breaks = data_legend$method_name, labels = data_legend$method_name) +
  facet_wrap(classification ~ genelist_size_label, ncol = 2) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE),
         shape = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = "RMSE", y = "Geneset size") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 9, margin = margin(2,1,2,1)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )
print(p)

ggsave(p, file = paste0(output_dir, "/null_simulations_bins__rmse.pdf"), width = 3.5, height = 6.5)


### print stats; average RMSE for each method
print(y |> group_by(method_name) |> summarise(mean_RMSE = mean(rmse_ngenes), .groups = "drop"))
