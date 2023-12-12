
The Geneset Ordinal Association Test (GOAT) is a parameter-free permutation-based algorithm for geneset enrichment analysis of preranked genelists. The full algorithm is computationally efficient and completes in the order of seconds and within 1 second when using precomputed null distributions. Validations using synthetic data show that estimated geneset p-values are well calibrated under the null hypothesis and invariant to geneset size. Application to various real-world proteomics and gene expression studies demonstrates that GOAT consistently identifies more significant Gene Ontology terms as compared to alternative methods.

The GOAT algorithm has not been published yet but a preprint is available, please cite it when using the early-access version of GOAT;

_Koopmans, F. (2023). GOAT: efficient and robust identification of geneset enrichment._ [https://doi.org/10.1101/2023.12.10.570979](https://doi.org/10.1101/2023.12.10.570979).


## Online tool

This GitHub repository contains the source-code for the GOAT R package and shows how to use it. Alternatively, check out the GOAT online tool at [https://ftwkoopmans.github.io/goat](https://ftwkoopmans.github.io/goat)


## Installation

The GOAT algoritm is implemented as an R package. Assuming R (at least version 4.1), RTools and RStudio are installed on your computer, execute the following R commands in a new RStudio session (i.e. close RStudio if currently opened, then start RStudio anew) to install GOAT;

```
install.packages("pak") # first, install the package manager that we'll use to install GOAT in the next step
pak::pkg_install("ftwkoopmans/goat", upgrade = FALSE) # latter parameter makes us skip optional updates

# optionally, only if you want to make use of fGSEA via the GOAT R package
pak::pkg_install("ctlab/fgsea", upgrade = FALSE)
```


## Quickstart

Below example assumes you've successfully installed the R package as per instructions above.

We'll use GOAT to identify enriched GO terms in one of the example proteomics datasets that is bundled with the GOAT R package. Besides the code comments that accompany each step in the example below, you can also check the respective R function documentation for further details (e.g. `?goat::load_genesets_go_bioconductor`). Below code will generate;

- an output table with geneset p-values
- a logfile that includes M&M text adapted to your settings
- basic visualizations of significant genesets

```
library(goat)
# download GO genesets  (only downloads once, caches data on your computer)
genesets_asis = load_genesets_go_bioconductor()
# example genelist, one of the datasets already bundled with the R package
# later examples will show how prepare your own data table / genelist
data("goat_example_datasets")
genelist = goat_example_datasets$`Hondius 2021:mass-spec:PMID33492460`

# importantly, we first filter all genesets to retain those relevant to your genelist
# for example, only retain genesets with at least 10 genes that are also in your genelist and remove those that contain more than 1500 genes in your genelist. Further, we remove genesets that contain more than half your genelist
# it is crucial that this function is applied for each genelist that you want to analyze because the resulting 'genesets_filtered' table is tailored to the specific genelist you are using at the moment
genesets_filtered = filter_genesets(genesets_asis, genelist, min_overlap = 10, max_overlap = 1500)

# apply the GOAT algorithm, then perform multiple testing correction using the Bonferroni method and consider proteins with adjusted p-value <= 0.05 as significant
# we here test for enrichment in the 'effectsize' column of your genelist. Alternatively, to test enrichment in the gene p-value dimension, replace 'effectsize' with 'pvalue'
# to use FDR for multiple testing correction instead of the stringent Bonferroni method, set padj_method = "BH"
result = test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)

# print the significant GO term counts, per geneset 'source' (CC/BP/MF), to console
print(result |> group_by(source) |> summarise(signif_count = sum(signif), .groups="drop"))

# store the results as an Excel table, and create a logfile that documents the GOAT settings you used
# this also includes M&M text that you can use in your publications.
#
# When running this example, change the output filename and path to some prexisting directory on your computer (use forward slashes in the file path)
save_genesets(result, genelist, filename = "C:/temp/goat.xlsx")

# generate lollipop charts for each GO domain (CC/BP/MF)
# When running this example, change the output directory to an existing directory on your computer
plot_lollipop(result, output_dir = "C:/temp/", topn = 50)
```


## Preparing your genelist as input for GOAT

**specification**

The expected format for your genelist is a data.frame (or 'tibble') that contains the following named columns;

- `gene` = required column with gene integer values (Human Entrez gene IDs is almost all cases)
- `symbol` = optional column with gene symbols that will be merged into the output tables
- `signif` = required column with logical (boolean) values that indicate which genes are foreground (TRUE) and background (FALSE). While not used in the GOAT algorithm, this is used by some general-purpose functions that e.g. prepare output tables
- `pvalue` required if you are using GOAT with `score_type='pvalue'` ; a column with finite numeric values
- `effectsize` required if you are using GOAT with `score_type='effectsize'` ; a column with finite numeric values


**example data**

First 3 lines of a table with only pvalue data that is used for GOAT (with `score_type='pvalue'`);

| gene  | symbol  | pvalue |
|-------|---------|--------|
| 348   |  APOE   |  0.01  |
| 335   |  APOA1  |  0.09  |
| 9948  |  WDR1   |  1     |

You may add the required `signif` column to indicate all genes with pvalue < 0.01 are considered foreground/significant using this R statement, assuming above genelist data.frame is called 'genelist'; 
`genelist$signif = genelist$pvalue < 0.01`


**but what if I only have gene symbols and no Entrez gene identifiers?**

The GOAT R package includes a convenience function to map gene symbols to human Entrez gene IDs; `symbol_to_entrez`

First, you need to download a data table from the www.genenames.org website; 

- download link: https://www.genenames.org/download/statistics-and-files/
- table: "Complete dataset download links" -->> "Complete HGNC approved dataset" -->> download the "TXT" table
- filename is typically something like hgnc_complete_set.txt

Next, you prepare a data.frame in R that holds your gene symbols and their respective pvalues and effectsizes, which might look like this;

| symbol | pvalue | effectsize |
|--------|--------|------------|
| APOE   |  0.01  |  2.7       |
| APOA1  |  0.09  |  -1.3      |
| WDR1   |  1     |  0.11      |

The following R code will map these genes to entrez and report the success/fail rate to console;

```
# update this file path to where you stored the HGNC data table
file_hgnc = "C:/data/hgnc_complete_set.txt"
hgnc = hgnc_lookuptable(file_hgnc)
genelist = symbol_to_entrez(genelist, hgnc)
# after ID mapping is done, we add a column named 'gene' to the genelist data.frame and populate it with the entrez gene IDs
genelist$gene = genelist$entrez_id 
```

As a next step you will need to remove failed gene mappings (i.e. no Entrez gene id was found) and remove redundant genes (same entrez_id) as shown in the next section.

**example R snippet to remove invalid rows**

Suppose that you prepared a data.frame named 'genelist' as in the above example and want to remove all rows that lack a valid value for 'gene', or remove duplicates (same gene ID on multiple rows). Then the following R snippet (that uses the dplyr package) is convenient;

```
genelist = genelist |>
  # remove rows that do not contain a numeric value for gene or pvalue
  filter(is.finite(gene) & is.finite(pvalue)) |>
  # sort the table by smallest/best pvalues on top
  arrange(pvalue) |>
  # retain only the first row for each unique gene
  distinct(gene, .keep_all = TRUE)
```


## Importing genesets

In typical use-cases, you can automatically download GO genesets via Bioconductor using the `load_genesets_go_bioconductor` function. Example @ quickstart section above.

Alternatively, several functions exist to import genesets from other sources. For each, please refer to the function documentation in R (e.g. issue the R command `?load_genesets_go_fromfile`) for download links to data files and example code.

- function `load_genesets_go_fromfile` : Load GO genesets directly from gene2go and .OBO files stored on your computer (e.g. to ensure you always use the exact same GO geneset definitions)

- function `load_genesets_syngo` : Load synaptic gene ontology (SynGO) genesets obtained from www.syngoportal.org

- function `load_genesets_gmtfile` : Load genesets from GMT files, a common file format for genesets (e.g. KEGG pathways obtained via MSigDB)


## Geneset simplification

We have implemented a basic geneset\*geneset similary measure that can be used to cluster genesets, allowing you to quickly identify redundant results (groups of genesets that have strong overlap).

Assuming you ran the quickstart example above, the following code will generate heatmap figures that may aid the interpretation of your GOAT results in case a large number of significant genesets were identified;

```
# this function generates geneset*geneset similarity matrices
clusters = cluster_genesets(result, genelist)
# find the subset of non-overlapping genesets. See the function documentation for tweaking these parameters using R command; ?reduce_genesets
# in brief, to collapse/reduce to a smaller list of genesets, lower simscore_threshold to e.g. 0.85 or 0.8  and/or  signifgenes_fraction to e.g. 0.8
result = reduce_genesets(clusters, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
# print signif geneset counts to console, before and after simplification
print(result |> filter(signif) |> count(source))
print(result |> filter(signif_and_reduced) |> count(source))

# generate heatmaps for each GO domain (CC/BP/MF). Again, don't forget to change the output directory to an existing directory on your computer and use forward slashes in the file path
plot_heatmap(clusters, "C:/temp/")
# repeat the lollipop plots made before, but now only for geneset that remain after simplification. See the function documentation for tweaking the plot, e.g. plotting only a subset of genesets/results
plot_lollipop(result, "C:/temp/", only_reduced = TRUE)
```
