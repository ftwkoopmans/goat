% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_genesets_go.R
\name{go_find_parents}
\alias{go_find_parents}
\title{from a directed edgelist (child-to-parent) to a recursive lookup of the final parent term / root}
\usage{
go_find_parents(child, parent)
}
\arguments{
\item{child}{array of GO term IDs that represent children}

\item{parent}{array of GO term IDs that represent respective parents}
}
\description{
from a directed edgelist (child-to-parent) to a recursive lookup of the final parent term / root
}
