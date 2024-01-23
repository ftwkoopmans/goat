
library(ggplot2)
devtools::load_all() # to reproduce this code on non-development system, replace with library(goat)
goat_print_version()


# hardcoded color-codings
clr = tibble::tibble(
  method = c("hypergeometric", "fisherexact_ease", "gsea_pvalue", "gsea_effectsize", "goat_pvalue", "goat_effectsize", "goat_pvalue_effectsize", "goat_effectsize_up", "goat_effectsize_down", "goat_effectsize_abs"),
  label = c("hypergeometric test", "Fisher-exact EASE", "GSEA p-value", "GSEA effectsize", "GOAT p-value", "GOAT effectsize", "GOAT p-value & effectsize", "GOAT effectsize up", "GOAT effectsize down", "GOAT effectsize abs"),
  clr = c("#69f0ae", "#ffd600", "#42a5f5", "#304ffe", "#ff7043", "#f50057", "#d50000", "#f50057", "#f50057", "#f50057")
)

# load datasets that are already bundled with the R package
data("goat_example_datasets")

# dataset_names_go_analysis = names(goat_example_datasets)
dataset_names_go_analysis = grep("26931375", names(goat_example_datasets), ignore.case = T, value = T, invert = T)

genesets_syngo = load_genesets_syngo("C:/DATA/SynGO_bulk_download_release_20231201/syngo_ontologies.xlsx", gene_database = "entrez")
# genesets_syngo = load_genesets_syngo("C:/DATA/SynGO_bulk_download_release_20210225/syngo_ontologies.xlsx", gene_database = "entrez")
genesets = load_genesets_go_fromfile(
  file_gene2go = "C:/VU/projects/Frank - GOAT (2022)/genesets/gene2go_2024-01-01.gz",
  file_goobo = "C:/VU/projects/Frank - GOAT (2022)/genesets/go_2024-01-01.obo"
)
# genesets = load_genesets_go_bioconductor() # the version of GO loaded here is determined by your Bioconductor installation

results_generic = NULL
padj_method = "bonferroni"

for(n in dataset_names_go_analysis) {
  genelist = goat_example_datasets[[n]]
  genesets_filtered = filter_genesets(genesets, genelist, min_overlap = 10L, max_overlap = 1500L, max_overlap_fraction = 0.5, min_signif = NA, max_size = NA, dedupe = FALSE)
  n_result = bind_rows(
    test_genesets(genesets_filtered, genelist, method = "goat", score_type = "pvalue", padj_method = padj_method, padj_cutoff = 0.05) |> mutate(method = "goat_pvalue"),
    test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = padj_method, padj_cutoff = 0.05) |> mutate(method = "goat_effectsize"),
    test_genesets(genesets_filtered, genelist, method = "hypergeometric", require_nsignif = 3L, padj_method = padj_method, padj_cutoff = 0.05) |> mutate(method = "hypergeometric"),      # note; prefiltering
    test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "pvalue", padj_method = padj_method, padj_cutoff = 0.05) |> mutate(method = "gsea_pvalue"),
    test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "effectsize", padj_method = padj_method, padj_cutoff = 0.05) |> mutate(method = "gsea_effectsize")
  )
  results_generic = bind_rows(results_generic, n_result |> mutate(dataset_label = n))


  pdf(paste0("C:/temp/goat_null-distribution_skew-normal-fit__", gsub("[^a-zA-Z0-9 _-]", " ", n), ".pdf"))
  hist(rankscore_fixed_order(nrow(genelist)), xlab = "Gene score", main = "")
  x = test_genesets(genesets_filtered, genelist, method = "goat_bootstrap", score_type = "pvalue", padj_method = padj_method, padj_cutoff = 0.05, verbose = TRUE)
  dev.off()
  Sys.sleep(0.5)
  graphics.off()
}




################################################################################################################################################################################
tib_plot = results_generic |>
  filter(dataset_label %in% dataset_names_go_analysis) |>
  mutate(
    dataset_label = sub(":", "\n", sub(":PMID\\d+", "", dataset_label, ignore.case = TRUE)),
    source = gsub("_", " ", source)
  ) |>
  group_by(dataset_label, source, method) |>
  summarise(n = sum(signif), .groups = "drop") |>
  mutate(method = factor(method, levels = rev(intersect(clr$method, unique(method))) ) )

clr_subset = clr |> filter(method %in% tib_plot$method)

p = ggplot(tib_plot, aes(x=n, y=source, fill = method, label = n)) +
  geom_col(position = ggplot2::position_dodge(.9)) +
  geom_text(position = ggplot2::position_dodge(.9), hjust = -0.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)), n.breaks = 3) +
  scale_fill_manual(values = tibble::deframe(clr_subset |> select(method, clr)), breaks = clr_subset$method, labels = clr_subset$label) +
  facet_wrap(.~dataset_label, nrow = 1, scales = "free_x") +
  labs(x = "significant genesets", y = "") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(0.75, "line"),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  )

print(p)
ggsave(p, file = "C:/temp/realworld_compare_signif_count.pdf", width = 8, height = 4.5)
Sys.sleep(1)
graphics.off()
# as SVG for website
svglite::svglite(filename = "C:/temp/barplot_signif_count.svg", width = 8, height = 4.5)


# total number of significant GO terms (across all CC/MF/BP domains) per dataset
tmp = results_generic |>
  filter(dataset_label %in% dataset_names_go_analysis) |>
  group_by(dataset_label, method) |>
  summarise(n = sum(signif), .groups = "drop")
# report the % increase of GOAT significant genesets as compared to GSEA  @  effectsize
tmp |>
  group_by(dataset_label) |>
  mutate(
    gain_gsea_n = n - n[method == "gsea_effectsize"],
    gain_gsea_pct = round((n - n[method == "gsea_effectsize"]) / n[method == "gsea_effectsize"] * 100, 1)
  ) |>
  ungroup() |>
  filter(method == "goat_effectsize") |>
  print(n = Inf)
# report the % increase of GOAT significant genesets as compared to GSEA  @  pvalues
tmp |>
  group_by(dataset_label) |>
  mutate(gain_gsea = round((n - n[method == "gsea_pvalue"]) / n[method == "gsea_pvalue"] * 100, 1) ) |>
  ungroup() |>
  filter(method == "goat_pvalue") |>
  print(n = Inf)
# report the % increase of GOAT significant genesets as compared to ORA  @  pvalues
tmp2 = tmp |>
  group_by(dataset_label) |>
  mutate(
    gain_ora_n = n - n[method == "hypergeometric"],
    gain_ora_pct = round((n - n[method == "hypergeometric"]) / n[method == "hypergeometric"] * 100, 1)
  ) |>
  ungroup()
tmp2 |> filter(method == "goat_pvalue") |> print(n = Inf)
tmp2 |> filter(method == "goat_effectsize") |> print(n = Inf)


lapply(goat_example_datasets, nrow)
lapply(goat_example_datasets, function(x)sum(x$signif))
lapply(goat_example_datasets, function(x)sum(x$signif)/nrow(x))



# for every dataset, plot top25 GOAT p-value datasets


## only for 2 * Colameo datasets
pdf("C:/temp/realworld_compare_scatterplot.pdf", 6, 15)
for(dataset_name in grep("Colameo 2021", unique(results_generic$dataset_label), value = T, ignore.case = T)) {

  tib_plot_filter = results_generic |>
    filter(dataset_label == dataset_name & method == "goat_pvalue") |>
    mutate(dataset_label = sub(":", ":\n", sub(":PMID\\d+", "", dataset_label, ignore.case = TRUE))) |>
    arrange(pvalue) |>
    group_by(source) |>
    mutate(myrank = 1:n()) |>
    filter(myrank <= 25) |>
    ungroup()

  tib_plot = results_generic |>
    filter(dataset_label == dataset_name & method %in% c("hypergeometric", "gsea_pvalue", "goat_pvalue") & id %in% tib_plot_filter$id) |>
    mutate(id = factor(id, levels = rev(tib_plot_filter$id)))

  clr_subset = clr |> filter(method %in% tib_plot$method)

  p = ggplot(tib_plot, aes(x = minlog10_fixzero(pvalue_adjust), y = id, size = ngenes, colour = method)) +
    geom_point(alpha = 0.5) +
    scale_fill_manual(values = tibble::deframe(clr_subset |> select(method, clr)), breaks = clr_subset$method, labels = clr_subset$label, aesthetics = c("fill", "colour")) +
    scale_y_discrete(breaks = levels(tib_plot$id), labels = string_trunc_right(tib_plot$name[match(levels(tib_plot$id), tib_plot$id)], width = 50)) +
    scale_size(range = c(2,5), trans = "log10", name = "gene count") +
    guides(colour = guide_legend(direction = "vertical", override.aes = list(alpha=1, size = 4)),
           size = guide_legend(direction = "vertical"), override.aes = list(alpha=1)) +
    facet_wrap( . ~ dataset_label + source, scales = "free", ncol = 1) +
    labs(x = "-log10(adjusted p-value)", y = "") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
  print(p)
}

dev.off()
Sys.sleep(1)
graphics.off()
################################################################################################################################################################################

writexl::write_xlsx(results_generic, "C:/temp/realworld_datasets_geneset_statistics.xlsx")
writexl::write_xlsx(results_generic |> filter(signif), "C:/temp/realworld_datasets_geneset_statistics_only_signif.xlsx")
save(results_generic, clr, file = "C:/temp/realworld_data.RData")

