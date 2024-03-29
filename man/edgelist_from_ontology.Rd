% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify_ontology.R
\name{edgelist_from_ontology}
\alias{edgelist_from_ontology}
\title{Construct DAG between genesets-of-interest as edgelist}
\usage{
edgelist_from_ontology(ids, genesets)
}
\arguments{
\item{ids}{shortlist of geneset IDs}

\item{genesets}{importantly, genesets should be the input genesets and not the "filtered" genesets because
only the former contains the complete ontological structure (parent/child links between genesets)}
}
\description{
this function reduces the complete DAG to a subset with only parameter genesets
}
