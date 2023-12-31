% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{parameters_prettyprint_length1}
\alias{parameters_prettyprint_length1}
\title{return a prettyprint string of all length 1 parameters that are string/numeric/logical/NA}
\usage{
parameters_prettyprint_length1(...)
}
\arguments{
\item{...}{arbitrary set of parameters}
}
\description{
return a prettyprint string of all length 1 parameters that are string/numeric/logical/NA
}
\examples{
\dontrun{
  parameters_prettyprint_length1(test1=1:2, test2=matrix(1:4,2,2), test3=data.frame(a=1),
                                 test4=c(a=1), test5=1, test6="a", test7=NA, test8=Inf)
}
}
