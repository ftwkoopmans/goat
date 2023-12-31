% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_genesets.R
\name{validate_genesets}
\alias{validate_genesets}
\title{validate that a geneset table is compatible with this R package}
\usage{
validate_genesets(genesets, require_signif = TRUE, check_gene_identiers = TRUE)
}
\arguments{
\item{genesets}{geneset tibble to validate, e.g. results from \code{load_genesets_go_bioconductor()}}

\item{require_signif}{should we check for "ngenes_signif" column ?}

\item{check_gene_identiers}{optionally skip the validation of gene identifiers}
}
\description{
validate that a geneset table is compatible with this R package
}
