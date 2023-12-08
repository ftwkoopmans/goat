

library(ggplot2)
devtools::load_all() # to reproduce this code on non-development system, replace with library(goat)
goat_print_version()

niter = 1000L

# load GO genesets and subset BP terms
genesets = load_genesets_go_fromfile(
  file_gene2go = "C:/VU/projects/Frank - GOAT (2022)/genesets/gene2go_2023-11-01.gz",
  file_goobo = "C:/VU/projects/Frank - GOAT (2022)/genesets/go_2023-11-01.obo"
) |>
  filter(source == "GO_BP")
# genesets = load_genesets_go_bioconductor() |> filter(source == "GO_BP") # the version of GO loaded here is determined by your Bioconductor installation

ugene = unique(unlist(genesets$genes))
# ensure our results are reproducible
set.seed(123)

result = NULL
for(iter_genelist_size in c(500L, 2000L, 10000L)) {
  for(iter_fraction_signif in c(0.2, 0.1, 0.01)) {
    for(i in 1L:niter) {
      # random genelist
      i_genelist = tibble::tibble(gene = sample(ugene, size = iter_genelist_size, replace = FALSE)) |>
        # define first N genes as significant (randomized genelist so ordering / which genes we define as signif doesn't matter)
        mutate(signif = seq_len(iter_genelist_size) <= iter_genelist_size * iter_fraction_signif)

      # we have to update the genesets table to apply filtering AND ensure the 'ngenes_signif' column is updated

      # ! non-standard setting for GOAT workflows, but specifically for this simulation; allow just 5 overlapping genes with dataset (e.g. otherwise possibly few terms @ genelist of size=500)
      i_genesets_filtered = filter_genesets(genesets, i_genelist, min_overlap = 5L, max_overlap = 1500L, min_signif = NA, max_size = NA, dedupe = FALSE)

      # now apply geneset filtering
      result = bind_rows(
        result,
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 0L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 0L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 1L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 1L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 2L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 2L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 3L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 3L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 4L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 4L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 5L)  |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 5L),
        test_genesets_hypergeometric(i_genesets_filtered, i_genelist, require_nsignif = 10L) |> select(id, ngenes, ngenes_signif, pvalue) |> filter(is.finite(pvalue)) |> mutate(iter = i, genelist_size = iter_genelist_size, fraction_signif = iter_fraction_signif, require_nsignif = 10L)
      )
    }

  }
}




# cover entire range of p-values in 200 bins
tib_plot_p100 = NULL
for(k in seq(0.005, 1, by = 0.005)) {
  tib_plot_p100 = bind_rows(
    tib_plot_p100,
    result |>
      group_by(genelist_size, fraction_signif, require_nsignif) |>
      summarise(
        frac_pval = sum(pvalue <= k) / n(),  # mean signif count over bootstrap iterations
        n = n(),
        k = k,
        .groups = "drop"
      )
  )
}

# zoom in on strongest p-values in 200 bins; 0.05
tib_plot_p005 = NULL
for(k in seq(0.00025, 0.05, by = 0.00025)) {
  tib_plot_p005 = bind_rows(
    tib_plot_p005,
    result |>
      group_by(genelist_size, fraction_signif, require_nsignif) |>
      summarise(
        frac_pval = sum(pvalue <= k) / n(),  # mean signif count over bootstrap iterations
        n = n(),
        k = k,
        .groups = "drop"
      )
  )
}




myplot = function(x, use_log10 = F) {
  x = x |>
    arrange(fraction_signif) |>
    mutate(
      fraction_signif = paste0("foreground ", round(fraction_signif * 100), "%"),
      fraction_signif = factor(fraction_signif, levels = unique(fraction_signif))
    ) |>
    arrange(genelist_size) |>
    mutate(
      genelist_size = paste0("genelist n=", genelist_size),
      genelist_size = factor(genelist_size, levels = unique(genelist_size))
    )

  if(use_log10) {
    x$k = -log10(x$k)
    x$frac_pval = -log10(x$frac_pval)
  }

  ggplot(x, aes(x = k, y = frac_pval, group = require_nsignif, colour = factor(require_nsignif))) +
    geom_line(size = 0.75) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
    facet_wrap(fraction_signif~genelist_size) +
    guides(colour = guide_legend(title = "#signif required", nrow = 1, byrow = T)) +
    labs(x = ifelse(use_log10, "Expected -log10 p-value", "Expected p-value"),
         y = ifelse(use_log10, "Observed -log10 p-value", "Observed p-value")) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      strip.text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
}


pqq = myplot(tib_plot_p100, TRUE)

p100 = myplot(tib_plot_p100) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


p005 = myplot(tib_plot_p005) +
  scale_x_continuous(limits = c(0, 0.06), breaks = c(0.01, 0.03, 0.05)) +
  scale_y_continuous(limits = c(0, 0.06), breaks = c(0.01, 0.03, 0.05))

ggsave(pqq, file = "C:/temp/bootstrap_ORA_minsignif_qqplot.pdf", width = 6, height = 6)
ggsave(p100, file = "C:/temp/bootstrap_ORA_minsignif_p100.pdf", width = 6, height = 6)
ggsave(p005, file = "C:/temp/bootstrap_ORA_minsignif_p005.pdf", width = 6, height = 6)


