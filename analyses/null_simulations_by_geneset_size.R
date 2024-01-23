
library(ggplot2)
devtools::load_all() # to reproduce this code on non-development system, replace with library(goat)
goat_print_version()
output_dir = "C:/temp"  # change the output folder where plots will be stored to whatever directory on your computer (note the use of forward slashes)



# hardcoded color-codings across plots/figures
clr = tibble::tibble(
  method = c("hypergeometric", "fisherexact_ease", "gsea_pvalue", "gsea_effectsize", "goat_pvalue", "goat_effectsize", "goat_pvalue_effectsize", "goat_effectsize_up", "goat_effectsize_down", "goat_effectsize_abs"),
  label = c("hypergeometric test", "Fisher-exact EASE", "GSEA p-value", "GSEA effectsize", "GOAT p-value", "GOAT effectsize", "GOAT p-value & effectsize", "GOAT effectsize up", "GOAT effectsize down", "GOAT effectsize abs"),
  clr = c("#69f0ae", "#ffd600", "#42a5f5", "#304ffe", "#ff7043", "#f50057", "#d50000", "#f50057", "#f50057", "#f50057")
)


benchmark_bin_results = function(x) {
  x |>
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
    ) |>
    # summary statistics, reducing the huge table to only information we will plot
    group_by(ngenes, bin) |>
    summarise(
      count_genesets = n(),
      frac_minlog10 = mean(frac_minlog10),
      pvalue_minlog10 = mean(pvalue_minlog10),
      .groups = "drop"
    )
}


# setup input data; GO genesets  &  real-world dataset that is to be tested against randomized genesets
genesets = load_genesets_go_fromfile(
  file_gene2go = "C:/VU/projects/Frank - GOAT (2022)/genesets/gene2go_2023-11-01.gz",
  file_goobo = "C:/VU/projects/Frank - GOAT (2022)/genesets/go_2023-11-01.obo"
)
# genesets = load_genesets_go_bioconductor() # the version of GO loaded here is determined by your Bioconductor installation

# setup input data; we'll sample from a real-world genelist
data(goat_example_datasets)
genelist = goat_example_datasets$`Sahadevan 2021:RNA-seq:PMID34021139`

# identifiers for the random genesets we'll generate. typically some vector of length 100k ~ 500k
# the length of this vector defines the number of random genesets we'll generate
idseq = 1L:200000L # 200k geneset permutations

# ensure our results are reproducible
set.seed(123)

tib_plot_all = NULL
# for(genelist_size in c(100L, 500L, 1000L, 5000L)) {
for(genelist_size in c(250L, 1000L, 4000L, 12000L)) {
  stopifnot(genelist_size < nrow(genelist))

  cat("genelist_size =", genelist_size, "\n")

  ### generate the mock data

  # subset N random genes from the genelist
  mock_genelist = genelist |> slice(sample(seq_len(nrow(genelist)), size = genelist_size))

  # generate many random genesets of size k
  genesets_mock = NULL
  # for(mock_size in c(10L, 20L, 30L, 40L, 50L, 100L, 200L, 1000L)) {
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
                             mutate(genelist_size = genelist_size, method_name = "GOAT effectsize"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_up")) |>
                             mutate(genelist_size = genelist_size, method_name = "GOAT effectsize_up"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_down")) |>
                             mutate(genelist_size = genelist_size, method_name = "GOAT effectsize_down"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "goat_precomputed", score_type = "effectsize_abs")) |>
                             mutate(genelist_size = genelist_size, method_name = "GOAT effectsize_abs"))
  # GSEA
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 1000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=1000"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 1000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA effectsize\nniter=1000"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 10000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=10000"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 10000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA effectsize\nniter=10000"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "pvalue", nPermSimple = 50000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA p-value\nniter=50000"))
  tib_plot_all = bind_rows(tib_plot_all, benchmark_bin_results(test_genesets(genesets = genesets_mock, genelist = mock_genelist, method = "gsea", score_type = "effectsize", nPermSimple = 50000L)) |>
                             mutate(genelist_size = genelist_size, method_name = "GSEA effectsize\nniter=50000"))

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



# Plot only geneset bins where we have at least 25 datapoints to prevent extreme values,
# that are observed only once or twice out of 200 thousand (!) random draws,
# from showing up in the plot and causing strong bias in interpretation
# Outliers not supported by N observations are not indicative of an algorithm's performance.
# (note that with different set.seed, those 1-out-of-200k extremes will vary while the rest should not due to averaging across bootstrap iterations)
myplot = function(x, min_datapoints = 25) {
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

  ggplot(y, aes(frac_minlog10, pvalue_minlog10, colour = ngenes)) +
    geom_line(size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", size = 0.4) +
    facet_wrap(method_label ~ genelist_size_label, ncol = n_distinct(y$method_label), dir = "v") +
    guides(colour = guide_legend(title = "geneset size", nrow = 1, reverse = T, override.aes = list(size = 1))) +
    scale_colour_discrete(direction = -1) +
    labs(x = "Expected -log10 p-value", y = "Observed -log10 p-value") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 9, margin = margin(2,1,2,1)),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}


if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save(tib_plot_all, file = paste0(output_dir, "/bootstrap_genesets.rda"), compress = TRUE)
# load("C:/VU/projects/Frank - GOAT (2022)/manuscript/figures_v2_post-bugfix/bootstrap_genesets.rda")

# Figure 2: GOAT pvalue, GOAT effectsize, GSEA 1kiter pvalue, GSEA 1kiter effectsize, GSEA 50kiter pvalue, GSEA 50kiter effectsize
plot_names = c("GOAT p-value", "GOAT effectsize", "GSEA p-value\nniter=1000", "GSEA effectsize\nniter=1000", "GSEA p-value\nniter=50000", "GSEA effectsize\nniter=50000")
p = myplot(tib_plot_all |> filter(method_name %in% plot_names) |> arrange(match(method_name, plot_names)) )
p = p + facet_grid(vars(genelist_size_label), vars(method_label))
print(p)
ggsave(p, file = paste0(output_dir, "/bootstrap_geneset_bins__figure2.pdf"), width = 6.5, height = 5)

# Figure S3: all GOAT score types
p = myplot(tib_plot_all |> filter(grepl("GOAT", method_name)) )
print(p)
ggsave(p, file = paste0(output_dir, "/bootstrap_geneset_bins__GOAT-all-scoretypes.pdf"), width = 10, height = 10)

# Figure S4: GSEA 1k/10k/50k iter compared
p = myplot(tib_plot_all |> filter(grepl("GSEA", method_name)) )
print(p)
ggsave(p, file = paste0(output_dir, "/bootstrap_geneset_bins__GSEA-compare-niter.pdf"), width = 10, height = 10)

# Figure for website
p = myplot(tib_plot_all |> filter(method_name == "GOAT p-value")) +
  facet_wrap(method_label ~ genelist_size_label, nrow = 1)
svglite::svglite(filename = "C:/temp/null_simulations.svg", width = 8, height = 3)
print(p)

# Figure with only GOAT p-value and effectsize
p = myplot(tib_plot_all |> filter(method_name == "GOAT p-value" | method_name == "GOAT effectsize")) +
  facet_wrap(method_label ~ genelist_size_label, nrow = 2)
print(p)
ggsave(p, file = paste0(output_dir, "/bootstrap_geneset_bins__GOAT-main-figure.pdf"), width = 5, height = 3.5)


