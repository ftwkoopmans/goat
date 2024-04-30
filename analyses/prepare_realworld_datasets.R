
### parse downloaded data into RData format
# to reproduce: download files to your computer and update respective file paths. References to papers/tables @ manuscript M&M

library(goat) # also loads dplyr

# TODO: setup input/output paths
dir_input = "C:/DATA/GOAT_datasets" # directory where all datasets are stored
dir_output = "C:/DATA/GOAT_datasets" # directory where the RData file should be stored
hgnc = hgnc_idmap_table("C:/DATA/hgnc/hgnc_complete_set__2024-01-01.txt") # e.g. previously downloaded from https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt


goat_manuscript_datasets = list()


dea = tibble::as_tibble(data.table::fread(paste0(dir_input, "/PMID34396684_proteomics.csv"), data.table = F))
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
goat_manuscript_datasets[["Colameo 2021:mass-spec:PMID34396684"]] = genelist



dea = tibble::as_tibble(data.table::fread(paste0(dir_input, "/PMID34396684_transcriptomics.csv"), data.table = F))
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
goat_manuscript_datasets[["Colameo 2021:RNA-seq:PMID34396684"]] = genelist



# note; no multiple testing correction in GO figure @ paper, they probably didn't use a background either
dea = readxl::read_excel(paste0(dir_input, "/PMID34021139_tableS8.xlsx"), sheet = 4)
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
goat_manuscript_datasets[["Sahadevan 2021:RNA-seq:PMID34021139"]] = genelist



dea = readxl::read_excel(paste0(dir_input, "/PMID33087358_tableS2_orig.xlsx"), skip = 6,
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
goat_manuscript_datasets[["Higginbotham 2020:mass-spec:PMID33087358"]] = genelist



dea = readxl::read_excel(paste0(dir_input, "/PMID32424284_tables4_orig.xlsx"), skip = 5)
genelist = dea |>
  select(Protein, log2fc = beta, sd = sd_beta, pvalue = p, pvalue_adjust = `BH adjusted p`) |>
  mutate(symbol = toupper(gsub("[ ;,|].*", "", Protein))) |>
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
goat_manuscript_datasets[["Wingo 2020:mass-spec:PMID32424284"]] = genelist



dea = readxl::read_excel(paste0(dir_input, "/PMID33492460_tables8_orig.xlsx"))
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
goat_manuscript_datasets[["Hondius 2021:mass-spec:PMID33492460"]] = genelist



########## Klaasen et al. IP experiment, data processed with MS-DAP
# dataset/download available via MS-DAP GitHub repository; https://github.com/ftwkoopmans/msdap/tree/master/examples/data

# load the MS-DAP R package
# library(msdap)
#
# dataset = import_dataset_metamorpheus(path = paste0(dir_input, "/dataset_Klaassen2018_pmid26931375", protein_qval_threshold = 0.05, collapse_peptide_by = "sequence_modified")
# dataset = import_fasta(dataset, files = c(paste0(dir_input, "/dataset_Klaassen2018_pmid26931375/UP000000589_10090.fasta"), paste0(dir_input, "/dataset_Klaassen2018_pmid26931375/UP000000589_10090_additional.fasta")))
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
# write.table(x, paste0(dir_input, "/pmid26931375_msdap.tsv"), quote = F, sep = "\t", row.names = F)
##########

dea = readr::read_tsv(paste0(dir_input, "/pmid26931375_msdap.tsv"))
genelist = dea |>
  filter(!is.na(entrez_id)) |>
  distinct(entrez_id, .keep_all = T) |>
  mutate(gene = entrez_id)
table(genelist$signif)

goat_manuscript_datasets[["Klaassen 2016:IP mass-spec:PMID26931375"]] = genelist





########## the list of datasets is now done

# save RData file to disk
save(goat_manuscript_datasets, file = paste0(dir_output, "/goat_manuscript_datasets.rda"), compress = "xz")

## unfortunately we cannot bundle a lot of data due to CRAN 5MB limit restrictions; so we only bundle the very small Klaassen et al. IP dataset (15kb compressed)
shisa6_apms_klaassen = goat_manuscript_datasets$`Klaassen 2016:IP mass-spec:PMID26931375` |> select(protein_id, log2fc, effectsize, pvalue, pvalue_adjust, signif, symbol, gene)
# uncomment line below to update R package (in dev environment)
# usethis::use_data(shisa6_apms_klaassen, overwrite = T, compress = "xz")


### alternatively, if we could fit a few hundred more kb, include the wingo et al. example
# usethis::use_data(goat_manuscript_datasets, overwrite = T, compress = "xz")
# wingo_proteomics_pmid32424284 = goat_manuscript_datasets$`Wingo 2020:mass-spec:PMID32424284` |> select(Protein, log2fc, effectsize, pvalue, pvalue_adjust, signif, symbol, gene)
# usethis::use_data(wingo_proteomics_pmid32424284, overwrite = T, compress = "xz")




