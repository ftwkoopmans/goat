
require(parallel)
library(ggplot2)
library(goat) # also loads dplyr

# TODO: setup input/output paths
output_dir = getwd() # directory where output files should be stored
file_gene2go = "gene2go_2024-01-01.gz" # full path to NCBI gene2go file, e.g. previously downloaded from https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
file_goobo = "go_2024-01-01.obo" # full path go GO OBO file, e.g. previously downloaded from http://current.geneontology.org/ontology/go.obo
source("test_genesets_idea.R") # load the script that implements the iDEA function; update to the full path


### parameters
# optionally disable iDEA, a method that takes a LONG time to compute
include_idea = TRUE
# fGSEA multiprocessing doesn't work on Windows. For other systems, set to max available cores minus one (cap at 24 max)
nthread = ifelse(tolower(.Platform$OS.type) == "windows", 1L, as.integer(min(24, max(1, parallel::detectCores() - 1))) )
cat(nthread, "thread(s) for fGSEA\n")


goat_print_version()








########################################################################################################################

# helper function that times how many seconds a parameter function takes to evaluate (discards the respective function's result)
mytime = function(f) {
  time_start = Sys.time()
  tmp = f()
  as.numeric(difftime(Sys.time(), time_start, units = "secs"))
}

# load GO genesets from gene2go file
genesets = load_genesets_go_fromfile(file_gene2go, file_goobo)

# for reproducibility of downstream random sampling
set.seed(123)

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
  cat("sim_geneset_count: ", sim_geneset_count, "\n")
  # we here use a subset of N genes and M genesets
  iter_genelist = mock_genelist |> utils::head(n = sim_geneset_count * 5)
  iter_genesets = filter_genesets(genesets, genelist = iter_genelist, min_overlap = 10, max_overlap = 1500L, max_overlap_fraction = 0.5) |>
    slice(sample(1L:n(), size = sim_geneset_count))

  tmp_time_idea = NA
  if(include_idea) {
    tmp_time_idea = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "test_genesets_idea", return_pvalue_louis = TRUE)})
  }

  # perform geneset enrichment tests and keep track of computation time
  result = bind_rows(result, tibble::tibble(
    geneset_count_target = sim_geneset_count,
    geneset_count = nrow(iter_genesets),
    genelist_size = nrow(iter_genelist),
    # fisherexact_ease            = mytime(function(){tmp = test_genesets(iter_genesets, iter_genelist, method = "fisherexact_ease")}),
    # hypergeometric              = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "hypergeometric")}),
    goat_pvalue                 = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat", score_type = "pvalue")}),
    goat_effectsize             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat", score_type = "effectsize")}),
    goat_fitfunction_pvalue     = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat_fitfunction", score_type = "pvalue")}),
    goat_fitfunction_effectsize = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "goat_fitfunction", score_type = "effectsize")}),
    gsea_pvalue_10k             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "pvalue", parallel_threads = nthread, nPermSimple = 10000L)}),
    gsea_effectsize_10k         = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "effectsize", parallel_threads = nthread, nPermSimple = 10000L)}),
    gsea_pvalue_50k             = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "pvalue", parallel_threads = nthread, nPermSimple = 50000L)}),
    gsea_effectsize_50k         = mytime(function(){test_genesets(iter_genesets, iter_genelist, method = "gsea", score_type = "effectsize", parallel_threads = nthread, nPermSimple = 50000L)}),
    idea                        = tmp_time_idea
  ))
}


# hardcoded color-codings
lbl = c(
  goat_pvalue                 = "GOAT p-value - precomputed",
  goat_effectsize             = "GOAT effect size - precomputed",
  goat_fitfunction_pvalue     = "GOAT p-value - full algorithm",
  goat_fitfunction_effectsize = "GOAT effect size - full algorithm",
  gsea_pvalue_10k             = "GSEA p-value - 10000 iterations",
  gsea_effectsize_10k         = "GSEA effect size - 10000 iterations",
  gsea_pvalue_50k             = "GSEA p-value - 50000 iterations",
  gsea_effectsize_50k         = "GSEA effect size - 50000 iterations",
  idea                        = "iDEA"#,
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
  gsea_effectsize_50k         = "#304ffe",
  idea                        = "#C99D8E"#,
  # hypergeometric              = "#69f0ae"
)


tib_plot = result |>
  # to long-format
  tidyr::pivot_longer(cols = setdiff(colnames(result), c("geneset_count_target", "geneset_count", "genelist_size")), names_to = "method", values_to = "time") |>
  # arrange factor levels according to the color-lookup-table
  mutate(
    label = lbl[match(method, names(lbl))],
    method = factor(method, levels = intersect(names(lbl), unique(method))),
    label = factor(label, levels = intersect(lbl, unique(label)))
  )


myplot = function(x) {
  stopifnot(x$label %in% names(lbl))
  clr_values = setNames(clr, lbl)
  clr_values = clr_values[names(clr_values) %in% x$label]
  ggplot(x, aes(geneset_count, time, colour = label, fill = label, shape = I(ifelse(grepl("effect.{0,1}size",label), 25, ifelse(grepl("p.{0,1}value",label), 24, 16))) )) +
    geom_point(size = 2.5) +
    geom_line(size = 0.75) +
    ylim(c(0, max(x$time))) +
    scale_colour_manual(values = clr_values, aesthetics = c("colour", "fill")) +
    guides(colour = guide_legend(ncol = 2, byrow = T, override.aes = list(alpha = 1))) +
    labs(x = "Geneset count", y = "Computation time (seconds)") +
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



save(result, tib_plot, lbl, clr, file = paste0(output_dir, "/computation_time.rda"), compress = "xz")
# finally, construct the ggplot
ggsave(myplot(tib_plot |> filter(method != "idea")), file = paste0(output_dir, "/computation_time.pdf"), width = 5.5, height = 5)
# analogous, for iDEA
if(include_idea) {
  x = tib_plot |> filter(method == "idea") |> mutate(time = time / 3600) # time on scale of hours, not seconds
  ggsave(myplot(x) + labs(y = "Computation time (hours)"), file = paste0(output_dir, "/computation_time_idea.pdf"), width = 5.5, height = 4)
}

