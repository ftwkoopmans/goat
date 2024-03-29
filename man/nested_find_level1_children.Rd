% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify_ontology.R
\name{nested_find_level1_children}
\alias{nested_find_level1_children}
\title{Traverse DAG and recursively split parent term into child terms until none are larger than N genes}
\usage{
nested_find_level1_children(obj, threshold, result = list())
}
\arguments{
\item{obj}{list representing a nested DAG}

\item{threshold}{stop if \code{obj$ngenes <= threshold}}

\item{result}{list of resulting elements (when calling this function, use default empty list)}
}
\description{
Traverse DAG and recursively split parent term into child terms until none are larger than N genes
}
