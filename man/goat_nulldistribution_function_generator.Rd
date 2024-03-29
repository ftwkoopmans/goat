% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/goat_nulldistributions.R
\name{goat_nulldistribution_function_generator}
\alias{goat_nulldistribution_function_generator}
\title{Wrapper function for goat_nulldistribution_function, typically used to generate the Rdata file with precomputed null distribution that is bundled with GOAT}
\usage{
goat_nulldistribution_function_generator(
  all_geneset_sizes,
  niter,
  dir_output,
  overwrite = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{all_geneset_sizes}{vector of geneset sizes; integer values larger than 10 and smaller than half the genelist length (gene_scores length)}

\item{niter}{see \code{goat_nulldistribution_function}}

\item{dir_output}{the path to an existing directory where temporary and result files will be stored}

\item{overwrite}{should existing files be overwritten? When setting this to FALSE, be mindful that you should manually clear cached files in between analyses !}

\item{verbose}{see \code{goat_nulldistribution_function}}
}
\description{
Wrapper function for goat_nulldistribution_function, typically used to generate the Rdata file with precomputed null distribution that is bundled with GOAT
}
\examples{
\dontrun{
# example 1; precompute parameters for a huge range of genesets, 500k permuations
goat_nulldistribution_function_generator(
  seq(from=100, to=20000, by=1), 500000, "C:/temp/", overwrite=TRUE, verbose=TRUE
)
# example 2; 2million permutations, a few very different genelist lengths; review QC plots
goat_nulldistribution_function_generator(
  c(100, 500, 1000, 5000, 15000), 2000000, "C:/temp/", overwrite=TRUE, verbose=TRUE
)
}
}
