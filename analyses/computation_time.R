
library(ggplot2)
devtools::load_all() # to reproduce this code on non-development system, replace with library(goat)
goat_print_version()


mytime = function(f) {
  time_start = Sys.time()
  tmp = f()
  as.numeric(difftime(Sys.time(), time_start, units = "secs"))
}


# for reproducibility of downstream random sampling
set.seed(123)

# load all GO BP ontology terms as genesets
genesets = load_genesets_go_fromfile(
  file_gene2go = "C:/VU/projects/Frank - GOAT (2022)/genesets/gene2go_2023-11-01.gz",
  file_goobo = "C:/VU/projects/Frank - GOAT (2022)/genesets/go_2023-11-01.obo"
)
# genesets = load_genesets_go_bioconductor() # the version of GO loaded here is determined by your Bioconductor installation

# generate mock genelist. values don't matter, all that we're evaluating here is computation time
mock_genelist = tibble::tibble(gene = unique(unlist(genesets$genes))) |>
  head(n = 20000) |>
  mutate(
    pvalue = runif(n()),
    log2fc = rnorm(n()),
    effectsize = log2fc,
    signif = pvalue <= 0.1
  )

# print number of signif genes in our mock genelist
print(table(mock_genelist$signif))

result = NULL
# try various geneset counts so we can check how methods scale against the number of genesets to test
for(sim_geneset_count in c(500, 1000, 2000, 4000, 6000)) {
  # repeat analyses a few times to ensure no method/algorithm is disadvantaged by unlucky RNG
  for(iter in 1:5) {
    # we here use a subset of N genes and M genesets
    iter_genelist = mock_genelist |> utils::head(n = sim_geneset_count * 5)
    iter_genesets = filter_genesets(genesets, genelist = iter_genelist, min_overlap = 10, max_overlap = 1500L, max_overlap_fraction = 0.5) |>
      slice(sample(1L:n(), size = sim_geneset_count))

    # perform geneset enrichment tests and keep track of computation time
    result = bind_rows(result, tibble::tibble(
      iteration = iter,
      geneset_count_target = sim_geneset_count,
      geneset_count = nrow(iter_genesets),
      genelist_size = nrow(iter_genelist),
      # fisherexact_ease            = mytime(function(){tmp = test_genesets(iter_genesets, iter_genelist, method = "fisherexact_ease")}),
      # hypergeometric              = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "hypergeometric")}),
      goat_pvalue                 = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat", score_type = "pvalue")}),
      goat_effectsize             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat", score_type = "effectsize")}),
      goat_fitfunction_pvalue     = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat_fitfunction", score_type = "pvalue")}),
      goat_fitfunction_effectsize = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat_fitfunction", score_type = "effectsize")}),
      gsea_pvalue_10k             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "pvalue", parallel_threads = 0L, nPermSimple = 10000L)}),
      gsea_effectsize_10k         = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "effectsize", parallel_threads = 0L, nPermSimple = 10000L)}),
      gsea_pvalue_50k             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "pvalue", parallel_threads = 0L, nPermSimple = 50000L)}),
      gsea_effectsize_50k         = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "effectsize", parallel_threads = 0L, nPermSimple = 50000L)})
    ))
  }
}



# hardcoded color-codings
lbl = c(
  goat_pvalue                 = "GOAT p-value - precomputed",
  goat_effectsize             = "GOAT effectsize - precomputed",
  goat_fitfunction_pvalue     = "GOAT p-value - full algorithm",
  goat_fitfunction_effectsize = "GOAT effectsize - full algorithm",
  gsea_pvalue_10k             = "GSEA p-value - 10000 iterations",
  gsea_effectsize_10k         = "GSEA effectsize - 10000 iterations",
  gsea_pvalue_50k             = "GSEA p-value - 50000 iterations",
  gsea_effectsize_50k         = "GSEA effectsize - 50000 iterations"#,
  # hypergeometric              = "hypergeometric test"
)

clr = c(
  goat_pvalue                 = "#ff7043",
  goat_effectsize             = "#f50057",
  goat_fitfunction_pvalue     = "#FF97E3",
  goat_fitfunction_effectsize = "#F069D3",
  gsea_pvalue_10k             = "#C788F8",
  gsea_effectsize_10k         = "#A83CFA",
  gsea_pvalue_50k             = "#42a5f5",
  gsea_effectsize_50k         = "#304ffe"#,
  # hypergeometric              = "#69f0ae"
)


tib_plot = result |>
  # to long-format, then collapse into mean values
  tidyr::pivot_longer(cols = setdiff(colnames(result), c("iteration", "geneset_count_target", "geneset_count", "genelist_size")), names_to = "method", values_to = "time") |>
  group_by(geneset_count, method) |>
  summarise(
    timesd = sd(time),
    time = mean(time),
    .groups = "drop"
  ) |>
  # arrange factor levels according to the color-lookup-table
  mutate(
    label = lbl[match(method, names(lbl))],
    method = factor(method, levels = intersect(names(lbl), unique(method))),
    label = factor(label, levels = intersect(lbl, unique(label)))
  )


decorate_plot = function(p) {
  p +
    scale_colour_manual(values = setNames(clr, lbl), aesthetics = c("colour", "fill")) +
    guides(colour = guide_legend(ncol = 2, byrow = T, override.aes = list(alpha = 1))) +
    labs(x = "geneset count", y = "computation time (seconds)") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.title.y = element_text(hjust = 1),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      plot.margin = margin(2, 70, 2, 30, unit = "pt")
    )
}

set.seed(123)
# finally, construct the ggplot
p = ggplot(tib_plot, aes(geneset_count, time, colour = label, fill = label, shape = I(c(24,25)[1+grepl("effectsize",label)])))
ggsave(decorate_plot(p + geom_point(aes(x = geneset_count  + runif(length(geneset_count), -90, 90)), size = 3, alpha = 0.7)),
       file = "C:/temp/computation_time_scatterplot.pdf", width = 5.5, height = 5)
# variant with lines
ggsave(decorate_plot(p + geom_point(size = 2) + geom_line(size = 0.75)),
       file = "C:/temp/computation_time_lineplot.pdf", width = 5.5, height = 5)


