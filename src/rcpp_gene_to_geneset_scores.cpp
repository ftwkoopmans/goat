// Copyright 2023 Frank Koopmans
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_gene_to_geneset_scores(int n_gs, Rcpp::IntegerVector gs_index, Rcpp::IntegerVector gs_geneindex, Rcpp::NumericMatrix gene_score)
{

  // basic input validation
  if ((gs_index.size() == 0) || (gs_index.size() != gs_geneindex.size()) || (gene_score.nrow() == 0) || (gene_score.ncol() == 0))
  {
    Rcpp::stop("invalid input data");
  }

  // optional bounds checking;
  // if(std::max_element(gs_index.begin(), gs_index.end()) > n_gs) { Rcpp::stop("gs_index > n_gs"); }
  // if(std::max_element(gs_geneindex.begin(), gs_geneindex.end()) > gene_score.nrow()) { Rcpp::stop("gs_geneindex > gene_score.nrow()"); }

  int N = gs_index.size();
  int n_scores = gene_score.ncol();
  Rcpp::NumericMatrix result(n_gs, n_scores);
  // std::fill( result.begin(), result.end(), 0 ); // zero is already the default value, we can skip this

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < n_scores; j++)
    {
      // indices were computed in R, we have to subtract 1 for C++ indices used here
      result(gs_index[i] - 1, j) += gene_score(gs_geneindex[i] - 1, j);
    }
  }

  return result;
}
