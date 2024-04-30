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


#include <random>
#include <cmath>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// in previous versions we here added an additional Rcpp requirement for C++14, but this seems to be no longer needed
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export(rcpp_geneset_null)]]
Rcpp::NumericVector rcpp_geneset_null(Rcpp::NumericVector gene_scores, Rcpp::IntegerVector geneset_sizes, int max_geneset_size, int niter)
{
  // note; we assume gene_scores have been rescaled upstream and are typically between 0 ~ 1000

  // basic input validation
  if ((niter < 10000) || (niter > 5000000))
  {
    Rcpp::stop("rcpp_geneset_null, failed precondition; 10000 <= niter <= 5000000");
  }
  if ((max_geneset_size < 5) || (max_geneset_size > 50000))
  {
    Rcpp::stop("rcpp_geneset_null, failed precondition; 5 <= max_geneset_size <= 50000");
  }
  if ((gene_scores.size() == 0) || (gene_scores.size() > 50000))
  {
    Rcpp::stop("rcpp_geneset_null, failed precondition; 0 < gene_scores.size() <= 50000");
  }
  if (max_geneset_size >= gene_scores.size())
  {
    Rcpp::stop("rcpp_geneset_null, failed precondition; max_geneset_size < gene_scores.size()");
  }
  if ((geneset_sizes.size() == 0) || (geneset_sizes.size() > 10000))
  {
    Rcpp::stop("rcpp_geneset_null, failed precondition; 0 < geneset_sizes.size() <= 10000");
  }

  std::vector<int> geneset_sizes__casted = Rcpp::as<std::vector<int>>(geneset_sizes); // type-cast up-front (minor speed improvement, maybe)
  int ngeneset_sizes = geneset_sizes.size();
  Rcpp::NumericVector result(ngeneset_sizes * niter); // following input validation checks, the size of this vector is at most 10k*1m = 1e+10

#ifdef _OPENMP
#pragma omp parallel
{
#endif
  // init openMP thread-local data
  std::vector<double> shuffled_gene_score = Rcpp::as<std::vector<double>>(gene_scores); // type-cast & copy source data into thread-local vector that will be shuffled
  std::vector<double> iter_geneset_scores(max_geneset_size);
  std::mt19937_64 myrng(123); // set seed to 123  (note that this is not the same RNG engine as we use in R, but that's fine)
  double cumsum = 0;

#ifdef _OPENMP
#pragma omp for schedule(static, 100)
#endif
  // openMP parallel loop
  for (int iter = 0; iter < niter; iter++)
  {
    // for reproducible pseudo-random bootstrapping, always set the same precomputed seed per iteration
    myrng.seed(iter + 1);
    // we do not need to reset the score vector (e.g. shuffled_gene_score = gene_scores__casted); just keep shuffling same thread-local vector -->> still reproducible WITHIN SAME SETTINGS
    std::shuffle(shuffled_gene_score.begin(), shuffled_gene_score.end(), myrng);

    // local results with scores for each geneset size (which equals index+1)
    // perhaps we can vectorize cumsum using e.g. std::accumulate or Armadillo. pseudocode; iter_geneset_scores = cumsum(head(shuffled_gene_score, max_geneset_size)) / (1:max_geneset_size)  // we can precompute/recycle the latter
    // but impact is probably low; compiler should be able to optimize this quite well
    cumsum = 0;
    for (int i = 0; i < max_geneset_size; i++)
    {
      cumsum += shuffled_gene_score[i];            // random gene score value, sampled from input gene score vector WITHOUT REPLACE
      iter_geneset_scores[i] = cumsum / (i + 1.0); // +1 to go from index to geneset size & implicit cast to double type
    }

    // only store the geneset sizes we are interested in (parameter geneset_sizes)
    // result array represent a matrix where row=geneset size, col=bootstrap iteration, values are stored "by row"
    // (i.e. first niter values are the first geneset size bin, second niter set of values is the second geneset size bin, etc.)
    for (int j = 0; j < ngeneset_sizes; j++)
    {
      // result vector index = current_row + current_column * number_of_columns
      // score = from geneset-size-index j to geneset size (geneset_sizes__casted[j]) -->> minus one to go from size to index -->> subset geneset score vector
      result[iter + j * niter] = iter_geneset_scores[geneset_sizes__casted[j] - 1];
    }
  }

#ifdef _OPENMP
}
#endif

// note that we don"t have to cleanup any C++ variables; we only use std::vector and primitive types

return result;
}



//// Welford running variance implementation copied from John Cook's blog
// https://www.johndcook.com/blog/standard_deviation/
class RunningStat
{
public:
  RunningStat() {
    Clear();
  }

  void Clear()
  {
    m_n = 0;
    m_oldM = m_newM = m_oldS = m_newS = 0.0;
  }

  void Push(double x)
  {
    m_n++;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (m_n == 1)
    {
      m_oldM = m_newM = x;
      m_oldS = 0.0;
    }
    else
    {
      m_newM = m_oldM + (x - m_oldM) / m_n;
      m_newS = m_oldS + (x - m_oldM) * (x - m_newM);

      // set up for next iteration
      m_oldM = m_newM;
      m_oldS = m_newS;
    }
  }

  int NumDataValues() const
  {
    return m_n;
  }

  double Mean() const
  {
    return (m_n > 0) ? m_newM : 0.0;
  }

  double Variance() const
  {
    return ((m_n > 1) ? m_newS / (m_n - 1) : 0.0);
  }

  double StandardDeviation() const
  {
    return sqrt(Variance());
  }

private:
  int m_n;
  double m_oldM, m_newM, m_oldS, m_newS;
};



//// wrapper function for bootstrapping + fit distributions
// [[Rcpp::export(rcpp_null_distributions)]]
Rcpp::List rcpp_null_distributions(Rcpp::NumericVector gene_scores, Rcpp::IntegerVector geneset_sizes, int max_geneset_size, int niter)
{
  // generate null distributions. also perform input validation
  Rcpp::NumericVector result = rcpp_geneset_null(gene_scores, geneset_sizes, max_geneset_size, niter);
  // init result vectors
  int ngeneset_sizes = geneset_sizes.size();
  Rcpp::NumericVector result_sigma(ngeneset_sizes);

#ifdef _OPENMP
#pragma omp parallel
{
#endif
  // init openMP thread-local data
  RunningStat rs; // use the Welford running variance implementation from John Cook's blog

#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
  // openMP parallel loop
  for (int i = 0; i < ngeneset_sizes; i++)
  {
    rs.Clear();
    int index_start = niter * i;
    // int index_end = index_start + niter;
    for (int j = 0; j < niter; j++)
    {
      rs.Push(result[index_start + j]);
    }

    result_sigma[i] = rs.StandardDeviation();
  }

#ifdef _OPENMP
}
#endif

return Rcpp::List::create(Rcpp::Named("null") = result, Rcpp::Named("sigma") = result_sigma);
}



// [[Rcpp::export(rcpp_dsnorm_logsum)]]
double rcpp_dsnorm_logsum(Rcpp::NumericVector x_min_mean, int N, double xi, double sd)
{
  constexpr double pi = 3.14159265358979323846;
  const double m1 = 2.0 / sqrt(2.0 * pi);
  const double mu = m1 * (xi - 1.0 / xi);
  const double g = 2.0 / (xi + 1.0 / xi);
  const double sigma = sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0 / (xi * xi)) + 2.0 * (m1 * m1) - 1.0);
  const double g_log = std::log(g);
  const double sigma_log = std::log(sigma);
  const double sd_log = std::log(sd);
  double z = 0;
  double Xi = 0;
  double sumdens = 0;

  for (int i = 0; i < N; i++)
  {
    z = (x_min_mean[i] / sd) * sigma + mu; // (x_min_mean[i] / sd) = standardize x (we already subtracted mean upstream)
    Xi = xi;
    if (z < 0)
    {
      Xi = 1.0 / xi;
    }
    double temp = fabs(z / Xi);
    // skipping sigma because we standardized upstream
    double z_Xi_dnorm_log = -(M_LN_SQRT_2PI + 0.5 * temp * temp); // M_LN_SQRT_2PI is defined @ cmath.h
    // sigma_log = sigma from actual skew-normal, sd_log = standard deviation passed to this function as parameter
    double i_density_log = g_log + z_Xi_dnorm_log + sigma_log - sd_log;
    // this shouldn't overflow because N is at most 10^6 in our case (i.e. sum of all log density values should fit in a double)
    sumdens += i_density_log;
  }

  return -1.0 * sumdens;
}
