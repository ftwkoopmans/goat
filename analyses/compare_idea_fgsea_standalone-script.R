# Comparison between iDEA and fGSEA on previously prepared genelist (RNAseq summary stats) and genesets (GO).
# This script is a minimal implementation without other dependences other than iDEA and fGSEA (and dplyr).
# Results show that iDEA yields many missing values and finds substantially fewer enriched genesets (multiple testing correction being equal).
#
# author: Frank Koopmans (frank.koopmans@vu.nl)
# date: 02-02-2024


### Preparing the input RData file for this script, using the GOAT R package
# library(goat)
# genelist_asis = goat::goat_example_datasets$`Sahadevan 2021:RNA-seq:PMID34021139`
# genesets_asis = goat::load_genesets_go_fromfile(
#   file_gene2go = "C:/VU/projects/Frank - GOAT (2022)/genesets/gene2go_2024-01-01.gz",
#   file_goobo = "C:/VU/projects/Frank - GOAT (2022)/genesets/go_2024-01-01.obo"
# )
# genesets = goat::filter_genesets(genesets_asis, genelist_asis, min_overlap = 10L, max_overlap = 1500L, max_overlap_fraction = 0.5) |> select(source, id, name, genes, ngenes)
# genelist = genelist_asis |> select(gene, log2fc, pvalue, pvalue_adjust)
# save(genesets, genelist, file = "C:/temp/sahadevan_pmid34021139.rda")



require(iDEA)
require(fgsea)
library(dplyr)

# hardcoded setting for our workstation with 24 cores
PARAM_NCORE = 23

### load the genelist and GO genesets from RData file
# we already prepared from the GO database the subset of genesets that have at least 10 genes overlapping with the input genelist
load("sahadevan_pmid34021139.rda")
head(genelist) # referred to as "summary" in iDEA tutorial
head(genesets) # referred to as "annotation" in iDEA tutorial


### prepare a genelist to iDEA format
# this code follows the iDEA tutorial, with minor adaption that we fix "infinite se_beta estimates" following pvalue=1 in input
# note; we use numeric gene IDs. To avoid type conversion / indexing shenanigans, we prefix integer gene IDs with a string
genelist_idea = data.frame(beta = genelist$log2fc, beta_var = NA, pvalue = genelist$pvalue, row.names = paste0("gene_", genelist$gene))
genelist_idea$zscore = stats::qnorm(genelist_idea$pvalue / 2.0, lower.tail=FALSE)
genelist_idea$se_beta = abs(genelist_idea$beta / genelist_idea$zscore) # approximate the standard error of beta (log2fc)
genelist_idea$se_beta[!is.finite(genelist_idea$se_beta)] = 10 # update/fix; default high value where beta_se is missing, yields beta_var = 100
genelist_idea$se_beta[genelist_idea$se_beta < 10^-16] = 10^-16 # guard against "zero standard error" (e.g. at log2fc=0; 0/x=0) to ensure downstream iDEA code doesn't suffer from division-by-zero issues
genelist_idea$beta_var = genelist_idea$se_beta^2


### prepare genesets to iDEA format
# from a data.frame with a list-column of gene IDs, to a gene*geneset identity matrix
genesets_to_matrix = function(gs, ugene) {
  # note; we use numeric gene IDs. To avoid type conversion / indexing shenanigans, we prefix integer gene IDs with a string
  m = matrix(0L, nrow = length(ugene), ncol = nrow(gs), dimnames = list(paste0("gene_", ugene), gs$id))
  for(i in 1:nrow(gs)) {
    m[,i] = as.integer(ugene %in% gs$genes[[i]])
  }
  return(m)
}

genesets_idea = as.data.frame(genesets_to_matrix(genesets, genelist$gene))
# preview the data we prepared
print(head(genelist_idea, n=10))
print(genesets_idea[1:10, 1:5])


### code QC to double-check our conversion from a list-of-genesets to a gene identity matrix for iDEA is correct
# validate that the rownames align between the tables we prepared for iDEA
stopifnot(rownames(genesets_idea) == rownames(genelist_idea))
# validate that the gene*geneset identity matrix that we constructed has the exact same gene-per-geneset-count as the input data in list format
tmp = colSums(genesets_idea)
stopifnot(genesets$ngenes == tmp[match(genesets$id, names(tmp))])
# validate 10 random genesets in detail; exact same gene identifiers
for(tmpid in sample(genesets$id, size = 10)) {
  stopifnot(paste0("gene_", genesets$genes[[match(tmpid, genesets$id)]]) %in% rownames(genesets_idea)[genesets_idea[,tmpid] == 1])
}


### following the iDEA tutorial, perform geneset enrichment testing
# as first parameter, we should select only the 2 required input columns and enforce their order (i.e. this is a hardcoded requirement in iDEA)
# we set `min_precent_annot=0` because we already filtered genesets upstream (i.e. default in tutorial is 0.0025, which amounts to genesets with 35 genes in a genelist of 14k genes; we do not want to skip small genesets)
idea = iDEA::CreateiDEAObject(genelist_idea[,c("beta", "beta_var")], genesets_idea, max_var_beta = 101, min_precent_annot = 0, num_core = PARAM_NCORE)
idea = iDEA::iDEA.fit(
  idea,
  fit_noGS = FALSE,
  init_beta = NULL,
  init_tau = c(-2,0.5),
  min_degene = 0, # defaults to 5 in the iDEA tutorial; "min_degene: the threshold for the number of detected DE genes in summary statistics. For some of extremely cases, the method does not work stably when the number of detected DE genes is 0"
  em_iter = 15,
  mcmc_iter = 1000,
  fit.tol = 1e-5,
  modelVariant = F,
  verbose = TRUE
)
# fix the the estimated geneset p-values as per default iDEA tutorial workflow
idea = iDEA::iDEA.louis(idea)


### map the iDEA estimated geneset p-values back to the input table
genesets$pvalue = idea@gsea$pvalue[match(genesets$id, idea@gsea$annot_id)]
genesets$pvalue_louis = idea@gsea$pvalue_louis[match(genesets$id, idea@gsea$annot_id)]
cat(sprintf("%d/%d genesets that were tested with iDEA received a p-value\n", sum(is.finite(genesets$pvalue)), nrow(genesets)))
cat("iDEA pvalue and pvalue_louis are finite value for the same genesets;\n")
print(table(is.finite(genesets$pvalue) == is.finite(genesets$pvalue_louis)))

# we've tested 3 GO sources (CC, BP, MF); apply stringent Bonferroni adjustment to each
genesets = genesets |>
  group_by(source) |>
  mutate(pvalue_adjust = stats::p.adjust(pvalue, method = "bonferroni")) |>
  ungroup()

# final result: print number of significant genesets per GO source
cat(sprintf("iDEA; %d/%d genesets are significant\n", sum(is.finite(genesets$pvalue_adjust) & genesets$pvalue_adjust <= 0.05), sum(is.finite(genesets$pvalue_adjust))))
print( genesets |> filter(pvalue_adjust <= 0.05) |> count(source) |> arrange(source) )



########## for reference, apply fGSEA with default settings
set.seed(123)
fgsea = fgsea::fgsea(
  pathways = stats::setNames(genesets$genes, genesets$id), # named geneset list
  stats = array(genelist$log2fc, dimnames = list(genelist$gene)), # named gene score array
  scoreType = "std",
  minSize = 1L, # disable filter, already done upstream
  maxSize = 10000L, # disable filter, already done upstream
  nPermSimple = 10000, # number of permutations (default is 1000)
  nproc = PARAM_NCORE
) |>
  select(id = pathway, pvalue = pval) |> # only need these 2 columns from fGSEA
  left_join(genesets |> select(id, source), by="id") # add the geneset "source" (GO CC/BP/MF)

# analogous to above code for adjusting geneset p-values
fgsea = fgsea |>
  group_by(source) |>
  mutate(pvalue_adjust = stats::p.adjust(pvalue, method = "bonferroni")) |>
  ungroup()
cat(sprintf("fGSEA; %d/%d genesets are significant\n", sum(is.finite(fgsea$pvalue_adjust) & fgsea$pvalue_adjust <= 0.05), sum(is.finite(fgsea$pvalue_adjust))))
print( fgsea |> filter(pvalue_adjust <= 0.05) |> count(source) |> arrange(source) )



########## write results to disk
save(idea, fgsea, file = "sahadevan_pmid34021139__idea_fgsea_results.rda", compress = "xz")



########## plot uncorrected p-values

pdf("sahadevan_pmid34021139__scatterplot_idea-vs-fgsea.pdf")
## scatterplot -log10 p-values (only shows data points where both methods yield a non-NA pvalue)
plot(-log10(fgsea$pvalue), -log10(genesets$pvalue[match(fgsea$id, genesets$id)]), xlab="-log10 pvalue fGSEA", ylab="-log10 pvalue iDEA")
abline(0, 1, col = 2)
# analogous for "louis" p-values
plot(-log10(fgsea$pvalue), -log10(genesets$pvalue_louis[match(fgsea$id, genesets$id)]), xlab="-log10 pvalue fGSEA", ylab="-log10 pvalue_louis iDEA")
abline(0, 1, col = 2)
dev.off()

## code QC: plot values directly from the iDEA result object (i.e. double-check there are no issues with ID mapping code)
# plot(-log10(fgsea$pvalue), -log10(idea@gsea$pvalue[match(fgsea$id, idea@gsea$annot_id)]), xlab="-log10 pvalue fGSEA", ylab="-log10 pvalue iDEA")
# abline(0, 1, col = 2)

