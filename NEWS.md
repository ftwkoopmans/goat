# *News*

# goat 1.0 (2024-04-28)

First public release.

# goat 1.1 (2025-01-19)

* save_genesets() will also store an RData file so users can easily load a previous analysis and create plots thereof.

* support gene sets from various organisms (besides human) via load_genesets_go_fromfile() and load_genesets_go_bioconductor(). Note that one has to provide an input gene list that contains Entrez gene identifiers of the same species.

* download_genesets_goatrepo() will automatically check the GitHub repository for the latest GO release if no specific date/version is provided.
