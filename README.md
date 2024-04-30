
The Gene set Ordinal Association Test (GOAT) is a parameter-free permutation-based algorithm for gene set enrichment analysis of preranked gene lists. The full algorithm is computationally efficient and completes in the order of seconds and within 1 second when using precomputed null distributions. Validations using synthetic data show that estimated gene set p-values are well calibrated under the null hypothesis and invariant to gene set size. Application to various real-world proteomics and gene expression studies demonstrates that GOAT consistently identifies more significant Gene Ontology terms as compared to alternative methods.

The GOAT algorithm has not been published yet but a preprint is available, please cite it when using the early-access version of GOAT;

_Koopmans, F. (2023). GOAT: efficient and robust identification of gene set enrichment._ [https://doi.org/10.1101/2023.12.10.570979](https://doi.org/10.1101/2023.12.10.570979).


# Contents
 - [Online tool](#online-tool) 
 - [Installing the R package](#installing-the-r-package) 
 - [Quickstart](#quickstart) 
 - [Tutorial: using the GOAT R package](#tutorial-using-the-goat-r-package) 
 - [Preparing your gene list as input for GOAT](#preparing-your-gene-list-as-input-for-goat) 
 - [Gene set simplification](#gene-set-simplification) 


# Online tool

This GitHub repository contains the source-code for the GOAT R package and shows how to use it. Alternatively, check out the user-friendly GOAT online tool at [https://ftwkoopmans.github.io/goat/](https://ftwkoopmans.github.io/goat/)


# Installing the R package

## 1) Preparing system dependencies

### Windows

- Install R version 4.1.0 (or higher) 
- Install RTools, available from https://cran.r-project.org/bin/windows/Rtools/
  - the RTools version you download/install from this link should match your current R version
  - if you don't know your currently installed R version (e.g. you already have R and RStudio but not RTools); start RStudio and check the version reported in the Console (bottom-left of the screen)
- Reboot your computer
- Install RStudio 

### MacOS

- For Apple computers, we suggest installing R from https://mac.r-project.org
  - if your current R version is below 4.3, you might have to upgrade by following these instructions to install R (installation process hasn't been tested with R 4.1 or 4.2 yet on MacOS)
  - download the installer for the current version listed under "R-4.X-branch"
    - note the difference between ARM (M1/M2/M3) and Intel-based
  - don't use the devel/unstable
- Install the 2 Mandatory tools described in https://mac.r-project.org/tools/   (Xcode and Fortran compiler)
  - following the official instructions linked here, the versions of R and these tools will be compatible and you should be able to install R packages from source
- Install RStudio 
- Start RStudio and install the rstudioapi package before proceeding with GOAT R package installation; `install.packages("rstudioapi")`

potential issues;

- `Directory of lock file does not exist: /Users/<username>/Library/Caches/org.R-project.R/R/pkgcache/pkg` during installation suggests you lack write permission to (some part of) this path. Solution; check that this path exists (create if missing) and give your current OS user write permission (at least to the last 4 directories in this path)
- `ld: library not found for -lgfortran` if this is among your error messages, the Fortran compiler is not installed (see second step in MacOS install instructions to resolve this) 
- `Failed to build source package <name>` --> either Xcode/Fortran is not installed (step 2 in above instructions), or their version doesn't line up with your current R version

### Linux

R version 4.1.0 (or higher) and toolchains for compilation are required. Suggested installation steps that include all system dependencies that we identified, for all (recursive) R packages that GOAT depends upon;

- Fedora; https://cran.r-project.org/bin/linux/fedora
  - `sudo dnf install R libcurl-devel libssl-devel libxml2-devel gmp-devel glpk-devel`
- Ubuntu; https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html
  - add R package repository to your /etc/apt/sources.list file (per instructions in above link)
  - `sudo apt-get update`
  - `sudo apt-get install r-base r-base-dev libcurl4-openssl-dev openssl-dev libxml2-dev libgmp-dev libglpk-dev libfontconfig1-dev libfreetype6-dev libgmp3-dev`
    - if your system reports "Unable to locate package openssl-dev", replace `openssl-dev` with `libssl-dev`
- (optional) install RStudio and then `install.packages("rstudioapi")`
- (if prompted, install packages into a personal library)
- install the R.utils package; start R (`R` on commandline), then `install.packages('R.utils')`

if some of the above dependencies in Fedora/Ubuntu are missing, you may encounter issues such as;

- `eval: make: not found` --> toolchain for compilation is missing, execute the above installation steps
- `Failed to build source package <name>` --> scroll up a bit to find the suggested solution. e.g. `Missing 2 system packages. You'll probably need to install them manually:` (followed by missing packages that need to be installed via dnf/apt)



## 2) Install the GOAT R package

The GOAT R package is available from GitHub (latest updates) and CRAN (most recent major release). Installation should be performed in a new RStudio session; close RStudio if currently opened, then start RStudio anew.

### install the latest version from GitHub

```
# if needed, install the package manager that we'll use to install GOAT in the next step
if (!requireNamespace("pak", quietly = TRUE)) 
  install.packages("pak")

# install latest goat version and all optional dependencies, but skip non-essential updates from other packages
pak::pkg_install("ftwkoopmans/goat", dependencies = TRUE, upgrade = FALSE)
```

### Alternatively, install the package from CRAN

```
# use the BiocManager tool to install optional goat dependencies from Bioconductor 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "GO.db", "org.Hs.eg.db"))

# install goat from CRAN
install.packages("goat")
```



## troubleshooting

- if `pak::pkg_install("ftwkoopmans/goat", upgrade = FALSE)` yields error like `! Failed to build source package goat` then (most likely) the toolchain for compilation is missing. Follow above instructions for Windows/MacOS/Linux
- if the installation of optional dependencies (e.g. Bioconductor packages) cause problems, try to set the `pak::pkg_install` parameter `dependencies = NA` to skip these



# Quickstart

Just show me the code! Download test data, apply GOAT, store result table.

```
library(goat)

# TODO: change the output dir to some folder on your computer
# e.g. "~/data/goat" on unix systems, "C:/data/goat" on Windows, 
# or set to getwd() to write output to the current working directory
output_dir = getwd()

# download an example gene list
datasets = download_goat_manuscript_data(output_dir)
genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`

# download GO gene sets
genesets_asis = download_genesets_goatrepo(output_dir)

# filter gene sets for sufficient overlap with the gene list
genesets_filtered = filter_genesets(genesets_asis, genelist)

# apply GOAT; score_type indicates the input gene list should be ranked by respective 
# "effectsize". This is the recommended setting. Alternatively,this can be set to "pvalue" 
# as well but will yield fewer genesets (c.f. GOAT manuscript, real-world benchmarks)
result = test_genesets(genesets_filtered, genelist, method = "goat",
  score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)

# print first 10 rows of the result table
print(result |> select(source, name, ngenes, pvalue_adjust) |> utils::head(n=10))

# store results. This function will create 2 files:
# - an output table with gene set p-values
# - a logfile that includes Methods text adapted to your settings
save_genesets(result, genelist, filename = paste0(output_dir, "/goat.xlsx"))
```



# Tutorial: using the GOAT R package

This tutorial assumes you have successfully installed the R package as per instructions above.

We'll use GOAT to identify enriched GO terms in an example proteomics dataset. Besides the comments that accompany each step in the example below, you can also check the respective R function documentation for further details (e.g. `?goat::load_genesets_go_fromfile`). Below code will generate;

- an output table with gene set p-values
- a logfile that includes Methods text adapted to your settings
- basic visualizations of significant gene sets

## load the R package and set output directory

```
library(goat)

# TODO: change the output dir to some folder on your computer
# e.g. "~/data/goat" on unix systems, "C:/data/goat" on Windows, 
# or set to getwd() to write output to the current working directory
output_dir = getwd()
```

## gene list

A preranked gene list is here defined as a table of gene identifiers and their respective effect sizes and/or p-values that indicate association with some experimental condition (e.g. summary statistics from an OMICS-based study). Below code downloads datasets that were analyzed in the GOAT manuscript and loads 1 as the gene list to be used in this tutorial. Documentation in the next section ("Preparing your gene list as input for GOAT") shows how to prepare your own data table / gene list.

```
datasets = download_goat_manuscript_data(output_dir)
genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`
```

## gene sets

A gene set can be any set of genes of interest; it is typically defined as a set of genes that are known members of the same biological pathway, localized to the same (sub)cellular compartment, co-expressed under certain conditions or associated with some disorder as defined in a gene set database such as GO. In the GOAT R package (and online tool) one tests for enrichment of top-ranked genes in the input gene list against each gene set from some collection/database. 

In the example below we'll use the GO database. To use another gene set database, refer to the section "Importing gene sets".

```
genesets_asis = download_genesets_goatrepo(output_dir)
```

## filter gene sets

Importantly, we first filter all gene sets to retain those relevant to your gene list.
In this example we use default (recommended) parameters to only retain gene sets with at least 10 genes that are also in your `genelist` and remove those that contain more than 1500 (or half the size of your gene list, whichever is smallest) genes that overlap with your gene list.

It is crucial that this function is applied for each gene list that you want to analyze because the resulting `genesets_filtered` table is only valid for the `genelist` you are using at the moment.

```
genesets_filtered = filter_genesets(genesets_asis, genelist, min_overlap = 10, max_overlap = 1500)
```

## gene set enrichment testing

Apply the GOAT algorithm, then perform multiple testing correction using the Bonferroni method and consider proteins with adjusted p-value <= 0.05 as significant. We here test for enrichment in the 'effectsize' column of your gene list. Alternatively, to test enrichment in the gene p-value dimension, replace 'effectsize' with 'pvalue'. But do keep in mind that the latter setting will yield fewer genesets (c.f. GOAT manuscript, real-world benchmarks).

To use FDR for multiple testing correction instead of the stringent Bonferroni method, set `padj_method = "BH"`

```
result = test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)

# print the significant GO term counts, per geneset 'source' (CC/BP/MF), to console
print(result |> group_by(source) |> summarise(signif_count = sum(signif), .groups="drop"))

# print the top 25 significant genesets (for simplicity, only select key data columns)
print(result |> filter(signif) |> select(source, name, ngenes, pvalue_adjust), n = 25)
```

## save results

Store the results as an Excel table, and create a logfile that documents the GOAT settings you used.
This logfile also includes Methods text tailored to your settings that is ready for use in scientific publications.

```
save_genesets(result, genelist, filename = paste0(output_dir, "goat.xlsx"))
```

## visualize results

Generate lollipop charts for each GO domain (CC/BP/MF), with gene set -log10 adjusted p-value on the x-axis and color-coding by gene set up/down-regulation.
Refer to the function documentation (`?plot_lollipop`) for alternative plot options (e.g. color-code by oddsratio on x-axis or create a barplot instead of a lollipop chart).

```
plot_lollipop(result, output_dir, topn = 50, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
```


# Preparing your gene list as input for GOAT

Gene identifiers used as input in the GOAT R package are Human NCBI Entrez gene IDs (other species currently not supported).

**specification**

The expected format for your gene list is a data.frame (or 'tibble') that contains the following named columns;

- `gene` = required column with gene integer values (Human Entrez gene IDs)
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

You may add the required `signif` column to indicate all genes with pvalue < 0.01 are considered foreground/significant using this R statement, assuming above gene list data.frame is called 'genelist'; 
`genelist$signif = genelist$pvalue < 0.01`


**but what if I only have gene symbols and no Entrez gene identifiers?**

The GOAT R package includes a convenience function to map gene symbols to human Entrez gene IDs; `symbol_to_entrez`

First, you need to download a data table from the www.genenames.org website; 

- download link: https://www.genenames.org/download/statistics-and-files/
- table: "Complete dataset download links" --> "Complete HGNC approved dataset" --> download the "TXT" table
- filename is typically something like hgnc_complete_set.txt

Next, you prepare a data.frame in R that holds your gene symbols and their respective pvalues and effectsizes, which might look like this;

| symbol | pvalue | effectsize |
|--------|--------|------------|
| APOE   |  0.01  |  2.7       |
| APOA1  |  0.09  |  -1.3      |
| WDR1   |  1     |  0.11      |

The following R code will map these genes to entrez and report the success/fail rate to console;

```
# TODO: update this file path to where you stored the HGNC data table (download link in above instructions)
file_hgnc = "C:/data/hgnc_complete_set.txt"
hgnc = hgnc_idmap_table(file_hgnc)
genelist = symbol_to_entrez(genelist, hgnc)
# the genelist table now contains Entrez gene IDs in the column 'gene'
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


# Importing gene sets

Several options are available to import gene sets from various sources. For each, please refer to the function documentation in R (e.g. issue the R command `?load_genesets_go_fromfile`) for further explanation and suggested download links to data files.

- function `download_genesets_goatrepo` : Automatically download GO gene sets from the GOAT GitHub repository. Data files are updated bianually.

- function `load_genesets_go_fromfile` : Load GO gene sets directly from gene2go and .OBO files stored on your computer. See the function help for download instructions.
This ensures you always use the exact same GO gene set definitions and gives you full control over the GO version used in your analyses.

- function `load_genesets_go_bioconductor` : Automatically download GO gene sets via Bioconductor. Note that the GO database version that is obtained may lag behind current data available from GO by years ! (depending on your currently installed R/Bioconductor version). The available GO data version is printed to console, keep an eye on this. If strongly outdated, either use the canonical `load_genesets_go_fromfile` function or update your R installation to the current release.

- function `load_genesets_syngo` : Load synaptic gene ontology (SynGO) gene set database obtained from www.syngoportal.org (see the function help for download instructions)

- function `load_genesets_gmtfile` : Load gene sets from GMT files, a common file format (e.g. KEGG pathways obtained via MSigDB). The [documentation section](https://ftwkoopmans.github.io/goat/docs#GMTgenesets) of the GOAT online webtool elaborates on how to obtain custom gene set databases (in GMT format).



# Gene set simplification

A basic geneset\*geneset similarity metric can be used to cluster gene sets, allowing you to quickly identify redundant results (groups of gene sets that have strong overlap).

Assuming you ran the quickstart example above, the following code will generate heatmap figures that may aid the interpretation of your GOAT results in case a large number of significant gene sets were identified;

```
## below code assumes the Quickstart example has been run (we here reuse variables 'result' and 'genelist')

# this function generates geneset*geneset similarity matrices
clusters = cluster_genesets(result, genelist)

# find the subset of non-overlapping gene sets. See the function documentation for tweaking these parameters using R command; ?reduce_genesets
# in brief, to collapse/reduce to a smaller list of genesets, lower simscore_threshold to e.g. 0.85 or 0.8  and/or  signifgenes_fraction to e.g. 0.8
result = reduce_genesets(clusters, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)

# print signif gene set counts to console, before and after simplification
print(result |> filter(signif) |> count(source))
print(result |> filter(signif_and_reduced) |> count(source))

# generate heatmaps for each GO domain (CC/BP/MF). Again, don't forget to change the output directory to an existing directory on your computer and use forward slashes in the file path
plot_heatmap(clusters, output_dir)

# repeat the lollipop plots made before, but now only for geneset that remain after simplification. See the function documentation for tweaking the plot, e.g. plotting only a subset of gene sets / results
plot_lollipop(result, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
```


Or try a treemap visualization of gene set testing results. Here, gene sets are grouped by their ontological structure, i.e. the parent/child relations between GO terms as described in the GO database.

```
## below code assumes the Quickstart example has been run (we here reuse variables 'result' and 'genesets_filtered')

# subset only significant CC terms
treemap_input = result |> filter(signif & source == "GO_CC")

# construct treemap data structures from gene set parent/child relations
tm = treemap_data(
  geneset_ids = treemap_input$id,
  genesets = genesets_filtered,
  genesets_test_result = treemap_input,
  ## simplify options;
  # "leaf_only" (most stringent, returns only leafs in the tree structure)
  # "prune_singletons" (remove parent terms that have exactly 1 child)
  # "pvalue" (remove parent terms where the child term p-value is at least 4 times better)
  # "none" (default; return all significant genesets that are not a "grouping term" in the treemap)
  simplify = "prune_singletons"
)

# plot the treemap
treemap_plot(tm$treemap_plotdata)
```

