# *News*

# goat 1.0 (2024-04-28)

* First public release.


# goat 1.1 (2025-01-19)

* save_genesets() will also store an RData file so users can easily load a previous analysis and create plots thereof.

* support gene sets from various organisms (besides human) via load_genesets_go_fromfile() and load_genesets_go_bioconductor(). Note that one has to provide an input gene list that contains Entrez gene identifiers of the same species.

* download_genesets_goatrepo() will automatically check the GitHub repository for the latest GO release if no specific date/version is provided.


# goat 1.1.2 (2025-02-22)

* bugfix: reduce_genesets()


# goat 1.1.3 (2025-08-27)

* if download functions fail they will yield a message instead of an error (compliance with CRAN policy). This affects functions; download_goat_manuscript_data() , download_genesets_goatrepo() , available_genesets_goatrepo()


# goat 1.1.4 (2026-01-24)

* plot_network() now returns the igraph object as an attribute
* reduce_genesets() bugfix: Add id to list of globals. Thanks @DavisVaughan !
