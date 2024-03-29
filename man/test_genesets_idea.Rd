% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_genesets_idea.R
\name{test_genesets_idea}
\alias{test_genesets_idea}
\title{Geneset testing with the iDEA R package (closely follows the iDEA tutorial)}
\usage{
test_genesets_idea(genesets, genelist, return_pvalue_louis)
}
\arguments{
\item{genesets}{filtered genesets}

\item{genelist}{user genelist, must contain columns "gene", "log2fc" and "pvalue"}

\item{return_pvalue_louis}{boolean, set \code{TRUE} to return "louis corrected" p-values from iDEA (analogous to iDEA tutorial)}
}
\value{
input \code{genesets} with additional columns for geneset p-value and score type. Because iDEA can return 2 types of geneset p-values,
multiple columns are included in this function's results; \code{pvalue_idea} = iDEA geneset p-value, \code{pvalue_louis_idea} = \code{pvalue_louis} from iDEA,
\code{pvalue} = respective copy of either former column, depending on user-parameter \code{return_pvalue_louis}
}
\description{
note that computation takes a long time with the iDEA method. On a high performance workstation computer with
12 cores / 24 threads it takes approximately 6 hours for a genelist of ~14000 genes to analyze all (filtered) GO
genesets. Using a regular PC this might take a multitude of hours.
}
\seealso{
\code{test_genesets}
}
