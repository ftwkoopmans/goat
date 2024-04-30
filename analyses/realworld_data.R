
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

padj_method = "bonferroni"
padj_cutoff = 0.05
# optionally disable iDEA, a method that takes a LONG time to compute
include_idea = TRUE




################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################

goat_print_version()

# load real-world OMICs datasets we prepared earlier from RData file
load(file_datasets)

# exclude the small and synapse-specific pulldown dataset from "typical OMICs datasets" for benchmarking
dataset_names_go_analysis = grep("26931375", names(goat_manuscript_datasets), ignore.case = T, value = T, invert = T)

# load GO genesets from gene2go file
genesets = load_genesets_go_fromfile(file_gene2go, file_goobo)

results = idea_timings = NULL
for(n in dataset_names_go_analysis) {
  file_cache = paste0(output_dir, "/cache__realworld__", gsub("[^a-zA-Z0-9 _-]", " ", n), ".rda")
  n_result = n_idea_timings = NULL
  if(file.exists(file_cache)) {
    load(file_cache)
  } else {
    cat("\nDATASET: ", n, "\n\n")
    gc(verbose = F, full = T)
    genelist = goat_manuscript_datasets[[n]]
    # filter genesets as per usual in GOAT workflow, removing genesets that have insufficient overlap with the input genelist
    genesets_filtered = filter_genesets(genesets, genelist, min_overlap = 10L, max_overlap = 1500L, max_overlap_fraction = 0.5, min_signif = NA, max_size = NA, dedupe = FALSE)
    # apply geneset enrichment algorithms, then remove columns that are not needed downstream to reduce RAM footprint of this script
    # (because we're applying many algorithms to many datasets; drop columns with all the gene constituents per geneset, etc.)
    n_result = bind_rows(
      test_genesets(genesets_filtered, genelist, method = "goat", score_type = "pvalue", padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "goat_pvalue"),
      test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "goat_effectsize"),
      test_genesets(genesets_filtered, genelist, method = "hypergeometric", require_nsignif = 3L, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "hypergeometric"), # note; ORA prefiltering
      test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "pvalue", parallel_threads = 1, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "gsea_pvalue"),
      test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "effectsize", parallel_threads = 1, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "gsea_effectsize"),
      # alternatively, also try GSEA with alternative weights. Ref; https://github.com/ctlab/fgsea/issues/45
      test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "pvalue", parallel_threads = 1, gseaParam = 2, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "gsea_pvalue_gseaParam2"),
      test_genesets(genesets_filtered, genelist, method = "gsea", score_type = "effectsize", parallel_threads = 1, gseaParam = 2, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif) |> mutate(method = "gsea_effectsize_gseaParam2")
    )

    if(include_idea) {
      # time one of the iDEA runs
      time_start = Sys.time()
      tmp_idea1 = test_genesets(genesets_filtered, genelist, method = "test_genesets_idea", return_pvalue_louis = TRUE, idea_variant = "default", verbose = TRUE, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif, ismissing_idea) |> mutate(method = "idea_default_louis")
      time_diff = as.numeric(difftime(Sys.time(), time_start, units = "secs"))
      # analogous for iDEA variants
      tmp_idea2 = test_genesets(genesets_filtered, genelist, method = "test_genesets_idea", return_pvalue_louis = TRUE, idea_variant = "rescale", verbose = TRUE, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif, ismissing_idea) |> mutate(method = "idea_rescale_louis")
      tmp_idea3 = test_genesets(genesets_filtered, genelist, method = "test_genesets_idea", return_pvalue_louis = TRUE, idea_variant = "beta_rescale", verbose = TRUE, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
        select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif, ismissing_idea) |> mutate(method = "idea_beta_rescale_louis")

      n_result = bind_rows(n_result, tmp_idea1, tmp_idea2, tmp_idea3)
      rm(tmp_idea1, tmp_idea2, tmp_idea3)
      n_idea_timings = tibble::tibble(dataset_label = n, genelist_length = nrow(genelist), geneset_count = nrow(genesets_filtered), time_sec = time_diff, time_hours = time_diff / 3600)

      ## debug;
      # x = test_genesets(genesets_filtered, genelist, method = "test_genesets_idea", return_pvalue_louis = TRUE, idea_variant = "beta_rescale_center", verbose = TRUE, padj_method = padj_method, padj_cutoff = padj_cutoff) |>
      #   select(id, source, ngenes_input, ngenes, ngenes_signif, pvalue, pvalue_adjust, signif, ismissing_idea) |> mutate(method = "idea_beta_rescale_center_louis")
    }

    # cache results
    save(n, n_result, n_idea_timings, file = file_cache, compress = "xz")

    # monitor progress; print |signif| per method (across all GO CC/MF/BP domains) to console
    print(n_result |> group_by(method) |> summarise(signif_count = sum(signif), .groups="drop"), n = Inf)

    # generate QC plots of the GOAT null distributions
    pdf(paste0(output_dir, "/realworld_goat_null-distribution_skew-normal-fit__", gsub("[^a-zA-Z0-9 _-]", " ", n), ".pdf"))
    hist(rankscore_fixed_order(nrow(genelist)), xlab = "Gene score", main = "")
    x = test_genesets(genesets_filtered, genelist, method = "goat_bootstrap", score_type = "pvalue", padj_method = padj_method, padj_cutoff = padj_cutoff, verbose = TRUE)
    dev.off()
    Sys.sleep(0.5)
    graphics.off()
  }

  # append dataset n to results
  results = bind_rows(results, n_result |> mutate(dataset_label = n))
  idea_timings = bind_rows(idea_timings, n_idea_timings)
}

if(include_idea) {
  cat("iDEA computation time per dataset:\n")
  print(idea_timings, n=Inf)
}

writexl::write_xlsx(results, paste0(output_dir, "/realworld_geneset_statistics.xlsx"))
save(results, clr, idea_timings, file = paste0(output_dir, "/realworld_data.rda"), compress = "xz")



################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################






########################################################################################################################
########################################## MAIN FIGURE: BARPLOT SIGNIF COUNTS ##########################################
########################################################################################################################

mybarplot = function(x) {
  tib_plot = x |>
    mutate(
      dataset_label = sub(":", "\n", sub(":PMID\\d+", "", dataset_label, ignore.case = TRUE)),
      source = gsub("_", " ", source)
    ) |>
    group_by(dataset_label, source, method) |>
    summarise(n = sum(signif), .groups = "drop") |>
    mutate(method = factor(method, levels = rev(intersect(clr$method, unique(method))) ) )

  clr_subset = clr |> filter(method %in% tib_plot$method)

  ggplot(tib_plot, aes(x=n, y=source, fill = method, label = n)) +
    geom_col(position = ggplot2::position_dodge(.9)) +
    geom_text(position = ggplot2::position_dodge(.9), hjust = -0.15, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25)), n.breaks = 4) +
    scale_fill_manual(values = tibble::deframe(clr_subset |> select(method, clr)), breaks = clr_subset$method, labels = clr_subset$label) +
    facet_wrap(.~dataset_label, nrow = 1, scales = "free_x") +
    labs(x = "significant genesets", y = "") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.75, "line"),
      # axis.ticks.x = element_blank(),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      panel.grid = element_blank()
    )
}

# ORA + multiple GSEA and iDEA variants
p = mybarplot(results)
ggsave(filename = paste0(output_dir, "/realworld_barplot_signif_count - everything.pdf"), plot = p, width = 8, height = 6)

# main figure; ORA + GOAT + updated/fixed iDEA variant  +  GSEA default/recommended
p = mybarplot(results |> filter(grepl("hypergeometric|goat|idea_default_louis|gsea_pvalue$|gsea_effectsize$", method)))
ggsave(filename = paste0(output_dir, "/realworld_barplot_signif_count - main figure.pdf"), plot = p, width = 8, height = 4.5)

p = mybarplot(results |> filter( ! grepl("hypergeometric|goat|idea_default_louis|gsea_pvalue$|gsea_effectsize$", method)))
ggsave(filename = paste0(output_dir, "/realworld_barplot_signif_count - suppl figure.pdf"), plot = p, width = 8, height = 4)




########################################################################################################################
############################################## PRINT SIGNIF INCREASE STATS #############################################
########################################################################################################################

x = results |> filter(grepl("idea", method)) |>
  group_by(dataset_label, method) |> summarise(n = n(), ismissing = sum(ismissing_idea %in% TRUE), pct = sprintf("%.1f", ismissing/n * 100), .groups = "drop") |>
  mutate(
    dataset_label = sub(":", "\n", sub(":PMID\\d+", "", dataset_label, ignore.case = TRUE)),
    method = factor(method, levels = c("idea_beta_rescale_louis", "idea_rescale_louis", "idea_default_louis")) # hardcode plot order
  ) |>
  arrange(dataset_label, method)
print(x)

p = ggplot(x, aes(y = method, x = ismissing / n * 100, label = pct)) +
  geom_col() +
  geom_text(hjust = -0.1) +
  facet_wrap(.~dataset_label) +
  coord_cartesian(c(0, 100)) +
  labs(x = "Fraction missing (%)", y = "") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank()
  )
ggsave(filename = paste0(output_dir, "/realworld_idea_missingness_barplot.pdf"), plot = p, width = 5.5, height = 4)




# total number of significant GO terms (across all CC/MF/BP domains) per dataset
tmp = results |>
  filter(dataset_label %in% dataset_names_go_analysis) |>
  group_by(dataset_label, method) |>
  summarise(n = sum(signif), .groups = "drop")
# report the % increase of GOAT significant genesets as compared to GSEA  @  effectsize
cat("*** gsea_effectsize\n")
tmp |>
  group_by(dataset_label) |>
  mutate(
    n_gsea = n[method == "gsea_effectsize"],
    gain_gsea_n = n - n[method == "gsea_effectsize"],
    gain_gsea_pct = round((n - n[method == "gsea_effectsize"]) / n[method == "gsea_effectsize"] * 100, 1)
  ) |>
  ungroup() |>
  filter(method == "goat_effectsize") |>
  print(n = Inf)
# report the % increase of GOAT significant genesets as compared to GSEA  @  pvalues
cat("*** gsea_pvalue\n")
tmp |>
  group_by(dataset_label) |>
  mutate(gain_gsea = round((n - n[method == "gsea_pvalue"]) / n[method == "gsea_pvalue"] * 100, 1) ) |>
  ungroup() |>
  filter(method == "goat_pvalue") |>
  print(n = Inf)
# report the % increase of GOAT significant genesets as compared to iDEA  @ effectsize
if(include_idea) {
  cat("*** idea_default_louis\n")
  tmp |>
    group_by(dataset_label) |>
    mutate(
      n_idea = n[method == "idea_default_louis"],
      gain_idea = n - n[method == "idea_default_louis"],
      gain_idea_pct = round((n - n[method == "idea_default_louis"]) / n[method == "idea_default_louis"] * 100, 1)
    ) |>    ungroup() |>
    filter(method == "goat_effectsize") |>
    print(n = Inf)
  cat("*** idea_rescale_louis\n")
  tmp |>
    group_by(dataset_label) |>
    mutate(
      n_idea = n[method == "idea_rescale_louis"],
      gain_idea = n - n[method == "idea_rescale_louis"],
      gain_idea_pct = round((n - n[method == "idea_rescale_louis"]) / n[method == "idea_rescale_louis"] * 100, 1)
    ) |>    ungroup() |>
    filter(method == "goat_effectsize") |>
    print(n = Inf)
  cat("*** idea_beta_rescale_louis\n")
  tmp |>
    group_by(dataset_label) |>
    mutate(
      n_idea_louis = n[method == "idea_beta_rescale_louis"],
      gain_idea_louis = n - n[method == "idea_beta_rescale_louis"],
      gain_idea_louis_pct = round((n - n[method == "idea_beta_rescale_louis"]) / n[method == "idea_beta_rescale_louis"] * 100, 1)
    ) |>
    ungroup() |>
    filter(method == "goat_effectsize") |>
    print(n = Inf)
}
# report the % increase of GOAT significant genesets as compared to ORA  @  pvalues
cat("*** hypergeometric\n")
tmp2 = tmp |>
  group_by(dataset_label) |>
  mutate(
    n_ora = n[method == "hypergeometric"],
    gain_ora_n = n - n[method == "hypergeometric"],
    gain_ora_pct = round((n - n[method == "hypergeometric"]) / n[method == "hypergeometric"] * 100, 1)
  ) |>
  ungroup()
tmp2 |> filter(method == "goat_pvalue") |> print(n = Inf)
tmp2 |> filter(method == "goat_effectsize") |> print(n = Inf)


lapply(goat_manuscript_datasets, nrow)
lapply(goat_manuscript_datasets, function(x)sum(x$signif))
lapply(goat_manuscript_datasets, function(x)sum(x$signif)/nrow(x))





########################################################################################################################
############################################# SCATTERPLOT GENESET P-VALUES #############################################
########################################################################################################################

# for every dataset, plot top25 GOAT p-value datasets
# only for 2 * Colameo datasets

pdf(paste0(output_dir, "/realworld_compare_topn_pvalmethods.pdf"), 6, 15)
for(dataset_name in grep("Colameo 2021", unique(results$dataset_label), value = T, ignore.case = T)) {

  tib_plot_filter = results |>
    filter(dataset_label == dataset_name & method == "goat_pvalue") |>
    mutate(dataset_label = sub(":", ":\n", sub(":PMID\\d+", "", dataset_label, ignore.case = TRUE))) |>
    arrange(pvalue) |>
    group_by(source) |>
    mutate(myrank = 1:n()) |>
    filter(myrank <= 25) |>
    ungroup()

  tib_plot = results |>
    filter(dataset_label == dataset_name & method %in% c("hypergeometric", "gsea_pvalue", "goat_pvalue") & id %in% tib_plot_filter$id) |>
    mutate(
      id = factor(id, levels = rev(tib_plot_filter$id)),
      name = genesets$name[match(id, genesets$id)]
    ) |>
    group_by(dataset_label, method) |>
    mutate(pvalue_adjust_minlog10 = minlog10_fixzero(pvalue_adjust, limit = NA)) |>
    ungroup()

  clr_subset = clr |> filter(method %in% tib_plot$method)

  p = ggplot(tib_plot, aes(x = pvalue_adjust_minlog10, y = id, size = ngenes, colour = method, shape = method)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(name = "method", values = tibble::deframe(clr_subset |> select(method, clr)) ) +
    scale_shape_discrete(name = "method") +
    scale_y_discrete(breaks = levels(tib_plot$id), labels = string_trunc_right(tib_plot$name[match(levels(tib_plot$id), tib_plot$id)], width = 50)) +
    scale_size(range = c(2,5), trans = "log10", name = "gene count") +
    guides(
      colour = guide_legend(direction = "vertical", ncol = 1, override.aes = list(size = 4)),
      shape = guide_legend(direction = "vertical", ncol = 1, override.aes = list(size = 4, ncol = 1)),
      size = guide_legend(direction = "vertical")) +
    facet_wrap( . ~ dataset_label + source, scales = "free", ncol = 1) +
    labs(x = "-log10(adjusted p-value)", y = "") +
    theme_bw() +
    theme(
      legend.position = "bottom"#,
      # legend.direction = "horizontal",
      # legend.box = "horizontal"
    )
  print(p)
}

dev.off()

















########################################################################################################################
################################################# INSPECT IDEA RESULTS #################################################
########################################################################################################################


# missed geneset p-values in iDEA vs GOAT; plot distributions (density lines) across datasets (facets/panels)
tmp = results |>
  # base data; geneset p-values from GOAT effectsize
  filter(method == "goat_effectsize") |>
  select(dataset_label, source, id, pvalue) |>
  # now add p-values from iDEA as separate column
  left_join(
    # remove "failed geneset p-value estimates" @ iDEA --> these will be NA in `tmp` table
    results |> filter(method == "idea_default_louis" & ! ismissing_idea %in% TRUE) |> select(dataset_label, source, id, pvalue_idea = pvalue),
    by = c("dataset_label", "source", "id")
  ) |>
  group_by(dataset_label) |>
  mutate(pvalue_minlog10 = minlog10_fixzero(pvalue, limit = NA)) |>
  ungroup()

# gather summary stat; total number of (filtered) genesets per dataset + count how many are missing in iDEA
tib_plot_stats = tmp |> group_by(dataset_label) |> summarise(ngeneset = n(), ngeneset_ideamiss = sum(is.na(pvalue_idea)), .groups = "drop")

# compose plot data
tib_plot = bind_rows(
  tmp |> mutate(classification = "all gene sets"),
  tmp |> filter(is.na(pvalue_idea)) |> mutate(classification = "missing in iDEA")
) |>
  mutate(classification = factor(classification, levels = c("all gene sets", "missing in iDEA"))) # enforce plot order


# helper function to update ggplot strip titles
create_strip_title = function(labels) {
  lapply(labels, function(lbl) {
    i = match(lbl, tib_plot_stats$dataset_label)
    paste0(
      # pretty-print label for dataset name (analogous to other plot)
      sub(":", " ", sub(":PMID\\d+", "", lbl, ignore.case = TRUE)),
      # add stat on # missing datasets in iDEA
      "\n", tib_plot_stats$ngeneset_ideamiss[i], "/", tib_plot_stats$ngeneset[i], " missing"
    )
  })
}


p = ggplot(tib_plot, aes(pvalue_minlog10, fill = classification, colour = classification)) +
  geom_density(alpha = 0.1) +
  scale_colour_manual(values = c("all genesets"="#455a64", "missing in iDEA"="#d50000"), aesthetics = c("colour", "fill")) +
  facet_wrap(.~dataset_label, ncol = 2, scales = "free", labeller = create_strip_title) +
  labs(x = "-log10 gene set p-value @ GOAT effect size", y = "Density") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(0.75, "line"),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  )
print(p)
ggsave(p, file = paste0(output_dir, "/realworld_characterize_idea_missing_genesets.pdf"), width = 6, height = 6)





################################################################################################################################################################################
# cross-compare existing methods; scatterplot fGSEA vs iDEA
for(plotspec in list(
  list(method="idea_default_louis", lbl="iDEA default"),
  list(method="idea_rescale_louis", lbl="iDEA rescaled"),
  list(method="idea_beta_rescale_louis", lbl="iDEA beta rescaled")
)) {
  tmp = results |>
    filter(method %in% c(plotspec$method, "gsea_effectsize")) |>
    select(id, pvalue, method, dataset_label) |>
    mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue))

  tib_plot = tmp |>
    tidyr::pivot_wider(id_cols = c("dataset_label", "id"), names_from = "method", values_from = "pvalue") |>
    # always name the iDEA pvalue column to "idea"
    rename(idea = !!sym(plotspec$method)) |>
    left_join(
      results |> filter(method == plotspec$method) |> mutate(in_idea = ! ismissing_idea %in% TRUE) |> select(dataset_label, id, in_idea),
      by = c("dataset_label", "id")
    ) |>
    mutate(missing_in_idea = factor( ! in_idea %in% TRUE, levels = c(FALSE, TRUE))) |>
    arrange(missing_in_idea)

  # compute R^2 using the subset of genesets that are in the top 25% of either method, while ignoring those with pvalue=1
  tib_plot_stats = tib_plot |>
    group_by(dataset_label) |>
    summarise(
      topn = n()*0.25, # number of genesets to use for current dataset; 25% "best"
      gsea_threshold = max(head(sort(gsea_effectsize, decreasing = F), n = topn)),
      idea_threshold = max(head(sort(idea, decreasing = F), n = topn)),
      nrow = sum(gsea_effectsize < 1 & is.finite(idea) & idea < 1 & (gsea_effectsize < gsea_threshold | idea < idea_threshold)),
      r2 = cor(minlog10_fixzero(gsea_effectsize[gsea_effectsize < 1 & is.finite(idea) & idea < 1 & (gsea_effectsize < gsea_threshold | idea < idea_threshold)], limit = NA),
               minlog10_fixzero(idea[           gsea_effectsize < 1 & is.finite(idea) & idea < 1 & (gsea_effectsize < gsea_threshold | idea < idea_threshold)], limit = NA),
               method = "pearson")^2,
      .groups = "drop"
    )
  # pretty-print label for the axis strip that includes the R^2
  mylabels = function(n) {
    x = sub(":", "\n", sub(":PMID\\d+", "", n, ignore.case = TRUE))
    sprintf("%s R^2=%.2f", x, tib_plot_stats$r2[match(n, tib_plot_stats$dataset_label)])
  }
  tib_plot$dataset_label = factor(tib_plot$dataset_label, levels = unique(tib_plot$dataset_label), labels = mylabels(unique(tib_plot$dataset_label)) )

  p = ggplot(tib_plot, aes(minlog10_fixzero(gsea_effectsize, limit = NA), minlog10_fixzero(idea, limit = NA), colour = missing_in_idea)) +
    geom_point(alpha = 0.66, shape = 16) +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
    facet_wrap(.~dataset_label, scales = "free") +
    guides(colour = guide_legend(title = paste("missing in", plotspec$lbl))) +
    scale_colour_manual(values = c("TRUE"="#d50000", "FALSE"="#455a64")) +
    labs(x = "-log10 p-value @ GSEA effect size", y = paste("-log10 p-value @", plotspec$lbl)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.75, "line"),
      axis.ticks.x = element_blank(),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 9),
      axis.title = element_text(size = 11)
    )

  print(p)
  ggsave(p, file = paste0(output_dir, "/realworld_scatterplot_gsea-vs-", plotspec$method, ".pdf"), width = 6, height = 6)
}


# analogous for GOAT vs GSEA (simplified variant of above plot, don't have to highlight missing)
tib_plot = results |>
  filter(method %in% c("goat_effectsize", "gsea_effectsize")) |>
  select(id, pvalue, method, dataset_label) |>
  tidyr::pivot_wider(id_cols = c("dataset_label", "id"), names_from = "method", values_from = "pvalue")
# compute R^2 using the subset of genesets that are in the top 25% of either method, while ignoring those with pvalue=1
tib_plot_stats = tib_plot |>
  group_by(dataset_label) |>
  summarise(
    topn = n()*0.25, # number of genesets to use for current dataset; 25% "best"
    gsea_threshold = max(head(sort(gsea_effectsize, decreasing = F), n = topn)),
    goat_threshold = max(head(sort(goat_effectsize, decreasing = F), n = topn)),
    r2 = cor(minlog10_fixzero(gsea_effectsize[gsea_effectsize < 1 & goat_effectsize < 1 & (gsea_effectsize < gsea_threshold | goat_effectsize < goat_threshold)], limit = NA),
             minlog10_fixzero(goat_effectsize[gsea_effectsize < 1 & goat_effectsize < 1 & (gsea_effectsize < gsea_threshold | goat_effectsize < goat_threshold)], limit = NA),
             method = "pearson")^2,
    .groups = "drop"
  )

# pretty-print label for the axis strip that includes the R^2
mylabels = function(n) {
  x = sub(":", "\n", sub(":PMID\\d+", "", n, ignore.case = TRUE))
  sprintf("%s R^2=%.2f", x, tib_plot_stats$r2[match(n, tib_plot_stats$dataset_label)])
}
tib_plot$dataset_label = factor(tib_plot$dataset_label, levels = unique(tib_plot$dataset_label), labels = mylabels(unique(tib_plot$dataset_label)) )

p = ggplot(tib_plot, aes(minlog10_fixzero(gsea_effectsize, limit = NA), minlog10_fixzero(goat_effectsize, limit = NA))) +
  geom_point(alpha = 0.66, shape = 16) +
  geom_abline(intercept = 0, slope = 1, colour = "black") +
  facet_wrap(.~dataset_label, scales = "free") +
  labs(x = "-log10 p-value @ GSEA effect size", y = "-log10 p-value @ GOAT effect size") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.75, "line"),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  )
print(p)
ggsave(p, file = paste0(output_dir, "/realworld_scatterplot_gsea-vs-goat.pdf"), width = 6, height = 6)



