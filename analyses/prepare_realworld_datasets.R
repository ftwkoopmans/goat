
### parse downloaded data into RData format

# to reproduce: download files to your computer and update respective file paths
# - see `hgnc_lookuptable` documentation for download instructions via genenames.org)

devtools::load_all() # to reproduce this code on non-development system, replace with library(goat)



hgnc = hgnc_lookuptable("C:/DATA/hgnc/hgnc_complete_set___2023-09-23.tsv")


goat_example_datasets = list()


dea = tibble::as_tibble(data.table::fread("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID34396684_proteomics.csv", data.table = F))
genelist = dea |>
  filter(Test == "Compartment") |>
  select(symbol = ID, log2fc = logFC, effectsize = t, pvalue = P.Value, pvalue_adjust = adj.P.Val) |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    signif = pvalue_adjust <= 10^-4
  ) |>
  filter(is.finite(gene) & is.finite(pvalue_adjust)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Colameo 2021:mass-spec:PMID34396684"]] = genelist



dea = tibble::as_tibble(data.table::fread("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID34396684_transcriptomics.csv", data.table = F))
genelist = dea |>
  filter(Test == "Compartment") |>
  select(symbol = names, log2fc = log2FC, effectsize = log2FC, pvalue = PValue, pvalue_adjust = FDR) |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    signif = pvalue_adjust <= 10^-4
  ) |>
  filter(is.finite(gene) & is.finite(pvalue_adjust)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Colameo 2021:RNA-seq:PMID34396684"]] = genelist



# note; no multiple testing correction in GO figure @ paper, they probably didn't use a background either
dea = readxl::read_excel("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID34021139_tableS8.xlsx", sheet = 4)
genelist = dea |>
  select(symbol = gene_name, log2fc = logFC, pvalue = PValue, pvalue_adjust = FDR) |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    effectsize = log2fc,
    signif = pvalue_adjust <= 0.05
  ) |>
  filter(is.finite(gene) & is.finite(pvalue_adjust)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Sahadevan 2021:RNA-seq:PMID34021139"]] = genelist



dea = readxl::read_excel("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID33087358_tableS2_orig.xlsx", skip = 6,
                         col_names = c("protein_id", "symbol", "module", "pvalue", "pvalue_adjust", "log2fc"))
genelist = dea |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    signif = pvalue_adjust <= 0.05,
    effectsize = log2fc
  ) |>
  filter(is.finite(pvalue) & is.finite(gene)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Higginbotham 2020:mass-spec:PMID33087358"]] = genelist



dea = readxl::read_excel("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID32424284_tables4_orig.xlsx", skip = 5)
genelist = dea |>
  select(Protein, log2fc = beta, effectsize = beta, pvalue = p, pvalue_adjust = `BH adjusted p`) |>
  mutate(symbol = toupper(gsub("[ ;,|].*", "", Protein))) |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    signif = pvalue_adjust <= 0.05
  ) |>
  filter(is.finite(pvalue) & is.finite(gene)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Wingo 2020:mass-spec:PMID32424284"]] = genelist



dea = readxl::read_excel("C:/VU/projects/Frank - GOAT (2022)/datasets/PMID33492460_tables8_orig.xlsx")
genelist = dea |>
  select(symbol = GENE, log2fc = `log2 FC Control vs Tangle`, effectsize = `log2 FC Control vs Tangle`, pvalue = `FDR Control vs Tangle`, pvalue_adjust = `FDR Control vs Tangle`) |>
  mutate(symbol = gsub(";.*", "", symbol)) |>
  symbol_to_entrez(hgnc = hgnc) |>
  mutate(
    symbol = hgnc_symbol,
    gene = entrez_id,
    signif = pvalue_adjust <= 0.05
  ) |>
  filter(is.finite(gene) & is.finite(pvalue_adjust)) |>
  arrange(pvalue) |>
  distinct(gene, .keep_all = T)
table(genelist$signif)
goat_example_datasets[["Hondius 2021:mass-spec:PMID33492460"]] = genelist



########## Klaasen et al. IP experiment, data processed with MS-DAP
# dataset/download available via MS-DAP GitHub repository; https://github.com/ftwkoopmans/msdap/tree/master/examples/data

# load the MS-DAP R package
# library(msdap)
#
# dataset = import_dataset_metamorpheus(path = "C:/VU/projects/Frank - GOAT (2022)/datasets/dataset_Klaassen2018_pmid26931375", protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified")
# dataset = import_fasta(dataset, files = c("C:/VU/projects/Frank - GOAT (2022)/datasets/dataset_Klaassen2018_pmid26931375/UP000000589_10090.fasta", "C:/VU/projects/Frank - GOAT (2022)/datasets/dataset_Klaassen2018_pmid26931375/UP000000589_10090_additional.fasta"))
# dataset = remove_proteins_by_name(dataset, regular_expression = "ig \\S+ chain|keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)") # optional
# dataset = sample_metadata_custom(dataset,
#                                  # a list of regular expressions to extract a short name and a group label from each full sample name
#                                  sample_property_regex = list(
#                                    shortname = c(".*CB_([A-Z]+.\\d).*", "\\1"), # regex dictates what part of the filename should be removed
#                                    group = c(".*CB_([A-Z]+)_.*", "\\1"),
#                                    fraction = c(".*\\D(\\d+)\\Dqtof.*", "\\1")
#                                  ),
#                                  # optionally, sort the groups in specified order
#                                  group_order = c("WT", "KO"))
# dataset = setup_contrasts(dataset, contrast_list = list(c("WT", "KO")))
#
#
# dataset = merge_fractionated_samples(dataset)
# dataset = dea(dataset, dea_algorithm = "deqms")
# dataset = differential_detect(dataset, min_peptides_observed = 2, min_samples_observed = 3, min_fraction_observed = 0.5, return_wide_format = FALSE)
# print_dataset_summary(dataset)
# plot_differential_detect(dataset)
#
# x = summarise_stats(dataset, return_dea = TRUE, return_diffdetect = TRUE, dea_logfc_as_effectsize = FALSE, diffdetect_zscore_threshold = 5, remove_ambiguous_proteingroups = FALSE)
#
# hgnc = hgnc_lookuptable(f = "C:/DATA/hgnc/hgnc_complete_set___2023-09-23.tsv")
# mgi = mgi_lookuptable(f = "C:/DATA/mgi/2023-07-15/MRK_SwissProt_TrEMBL.rpt")
# x = protein2gene_orthologs(x, hgnc, mgi)
# # print proteingroups that we cannot map to human genes
# print(x |> filter(is.na(hgnc_id)))
# write.table(x, "C:/VU/projects/Frank - GOAT (2022)/datasets/pmid26931375_msdap.tsv", quote = F, sep = "\t", row.names = F)
##########

dea = readr::read_tsv("C:/VU/projects/Frank - GOAT (2022)/datasets/pmid26931375_msdap.tsv")
genelist = dea |>
  filter(!is.na(entrez_id)) |>
  distinct(entrez_id, .keep_all = T) |>
  mutate(gene = entrez_id)
table(genelist$signif)

goat_example_datasets[["Klaassen 2016:IP mass-spec:PMID26931375"]] = genelist





# to reproduce: update file path
# save(goat_example_datasets, file = "C:/temp/goat_example_datasets.rda", compress = "xz")
# usethis::use_data(goat_example_datasets, overwrite = T, compress = "xz")
