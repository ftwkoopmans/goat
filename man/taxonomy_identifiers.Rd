% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{taxonomy_identifiers}
\alias{taxonomy_identifiers}
\title{Lookup table for taxonomy identifiers, their names and respective Bioconductor packages}
\usage{
taxonomy_identifiers()
}
\description{
To prevent adding a huge number of dependencies to this R package,
as the respective bioconductor_package also has to be added to the DESCRIPTION > Suggests,
we've only included the limited set of commonly used organisms.
}
