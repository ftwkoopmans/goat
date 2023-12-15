
#' compute rank^2 scores and rescale these between 0~1000 (with further precision captured by decimals)
#'
#' rank^2 can yield huge values, e.g. a large genelist of N=50000 genes would imply max gene score = 50000^2 = 2.5e+09.
#' Thus we rescale these scores between 0~1000 so downstream applications (e.g. sum of 10000 gene scores) don't explode to huge numbers.
#'
#' @param x data.frame to sort (i.e. the genelist)
#' @param sort1 numeric vector of length `nrow(x)` to sort by first/primarily (descending order, higher value gets better result score)
#' @param sort2 numeric vector of length `nrow(x)` to sort by for breaking ties (descending order, higher value gets better result score)
#' @param sort3 numeric vector of length `nrow(x)` to sort by for breaking ties (descending order, higher value gets better result score)
#' @param colname name for result column in `x`  (overwritten if already exists)
#' @examples \dontrun{
#'   x = data.frame(gene=c(1,2,3,4,5), pvalue=c(0.01, 1, 1, 0.1, 1),
#'                  effectsize=c(-2, 0.25, 0.5, 1, 0.25))
#'   print(x, row.names = FALSE)
#'   print(goat:::rankscore(x, sort1 = -1*x$pvalue, sort2 = abs(x$effectsize),
#'         sort3 = x$gene, colname="score") |>
#'     arrange(desc(score)), row.names = FALSE)
#' }
#' @return data.frame `x` with added column `colname`, containing gene scores between 0 and 1000
rankscore = function(x, sort1, sort2, sort3, colname) {
  i = NULL
  if(is.null(sort3)) {
    i = order(sort1, sort2, decreasing = TRUE)
  } else {
    i = order(sort1, sort2, sort3, decreasing = TRUE)
  }
  # add scores
  x[i,colname] = rev(rankscore_fixed_order(nrow(x)))
  return(x)
}



#' mock data on same scale as `rankscore`
#'
#' @param n genelist length
rankscore_fixed_order = function(n) {
  (1:n)^2 / (n^2 / 1000) # compute scaling factor prior to application to rank^2 scores
}



#' -log10 transform a vector of p-values, replacing zeros with some limit/threshold
#'
#' @examples \dontrun{
#'   pval = c(0, 10^-6, 0.001, 0.01, 1, NA, -Inf, Inf, NaN)
#'   cbind(
#'     input = pval,
#'     # default; replace zeros with typical R machine precision for doubles
#'     minlog10_default = minlog10_fixzero(pval),
#'     # alternatively, replace zero with lowest non-zero pvalue in input
#'     minlog10_limit_from_data minlog10_fixzero(pval, limit = NA)
#'   )
#' }
#' @param x p-value vector to transform to -log10
#' @param limit value to replace zero's in `x` with. Set NA to replace zero's in `x` with the smallest finite value in `x` (if there is none, defaults to 2.22e-16)
#' @export
minlog10_fixzero = function(x, limit = 2.22e-16) {
  stopifnot(length(limit) == 1 && (is.na(limit) || (is.numeric(limit) & is.finite(limit) & limit > 0)))
  stopifnot(is.na(x) | is.numeric(x)) # accept numbers, NA, NaN, Inf, -Inf. But not strings, booleans, factors, etc.

  # zero length input; return empty numeric
  if(length(x) == 0) return(numeric())

  # limit; if NA, use lowest non-zero value in input. But if no non-zero in input
  x_nonzero = x[is.finite(x) & x > 0]
  limit = ifelse(!is.finite(limit) && length(x_nonzero) > 0, min(x_nonzero), limit)
  limit = ifelse(!is.finite(limit), 2.22e-16, limit) # fallback for `x=0;limit=NA`; no x_nonzero to estimate from

  x[is.finite(x) & x == 0] = limit
  x[!is.finite(x)] = NA # i.e. replace NA, NaN, Inf, -Inf to NA prior to log transformation
  -log10(x)
}



#' Adjust p-values for all genesets, grouped by 'source' then adjust for the number of 'sources'
#'
#' @param genesets tibble with genesets, must contain column 'pvalue'
#' @param method method for multiple testing correction, must be any of `stats::p.adjust.methods`, e.g. "BH" or "bonferroni"
#' @param cutoff numeric cutoff value for adjusted p-value, `signif` column is set to TRUE for all values lesser-equals
#' @param correct_sources apply Bonferroni adjustment to all p-values according to the number of geneset sources that were tested. Boolean parameter, set TRUE to enable (default) or FALSE to disable
#' @export
padjust_genesets = function(genesets, method = "BH", cutoff = 0.01, correct_sources = TRUE) {
  pvalue = pvalue_adjust = NULL # fix invisible bindings R package NOTE
  stopifnot(length(genesets) > 0 && is.data.frame(genesets) && all(c("source", "pvalue") %in% colnames(genesets)))
  stopifnot(length(method) == 1 && method %in% stats::p.adjust.methods)
  stopifnot(length(cutoff) == 1 && is.finite(cutoff))
  stopifnot(length(correct_sources) == 1 && correct_sources %in% c(TRUE, FALSE))

  genesets = genesets |>
    group_by(source) |> # new group_by() will override preexisting groups (as here intended)
    mutate(pvalue_adjust = stats::p.adjust(pvalue, method = method)) |>
    ungroup()

  # correction for the number of sources that were evaluated
  if(correct_sources) {
    N = n_distinct(genesets$source)
    genesets = genesets |> mutate(pvalue_adjust = pvalue_adjust * N)
  }

  genesets |>
    mutate(
      pvalue_adjust = pmin(pvalue_adjust, 1), # cap p-value at max 1
      signif = is.finite(pvalue_adjust) & pvalue_adjust <= cutoff
    )
}



#' Density function for skew normal distribution
#'
#' @param x numeric vector of values
#' @param mean location parameter
#' @param sd scale parameter
#' @param xi skewness parameter
dsnorm = function(x, mean = 0, sd = 1, xi = 1) {
  dsnorm__std = function(x, xi) {
    if(length(xi) == 1) {
      xi = rep(xi, length(x))
    }
    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    z = x * sigma + mu
    g = 2 / (xi + 1 / xi)
    flag_up = z >= 0 # indicates upper half of the distribution
    XI = xi # flip xi for lower half
    XI[!flag_up] = 1 / XI[!flag_up]
    g * stats::dnorm(x = z/XI) * sigma
  }

  return(dsnorm__std(x = (x - mean)/sd, xi = xi) / sd)
}



#' Probability function for skew normal distribution
#'
#' @param q numeric vector of quantiles
#' @param mean location parameter
#' @param sd scale parameter
#' @param xi skewness parameter
psnorm_upper_tail = function(q, mean = 0, sd = 1, xi = 1) {
  # vectorized over xi
  psnorm_upper_tail__std = function(q, xi) {
    if(length(xi) == 1) {
      xi = rep(xi, length(q))
    }

    m1 = 2/sqrt(2 * pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
    z = q * sigma + mu
    g = 2 / (xi + 1 / xi)
    flag_up = z >= 0 # indicates upper half of the distribution
    XI = xi # flip xi for lower half
    XI[!flag_up] = 1 / XI[!flag_up]
    p = z # init as numeric vector of expected length
    # high precision for small pvalues
    p[flag_up] = g[flag_up] * XI[flag_up] * stats::pnorm(q = z[flag_up]/XI[flag_up], lower.tail = F)
    # precision for large p-values is capped near 1 (approx 1 - 10^-16), but doesn't matter
    p[!flag_up] = 1 - g[!flag_up] * XI[!flag_up] * stats::pnorm(q = z[!flag_up]/XI[!flag_up], lower.tail = T)

    p[p > 1] = 1
    return(p)
  }

  return( psnorm_upper_tail__std(q = (q - mean)/sd, xi = xi) )

  ### QC code
  # # validate precision; compare to pnorm @ xi = 1
  # q = seq(-12, 12, by = 0.1)
  # stopifnot( all.equal(stats::pnorm(q, lower.tail = F), psnorm_upper_tail__std(q, xi = 1), tolerance = 10^-14) )
  #
  # # validate correctness; compare to preexisting implementation
  # # here only using a limited range for q because reference lacks precision for large q (ours should be fine)
  # q = seq(-8, 8, by = 0.1)
  # for(xi in c(0.5, 0.8, 1, 1.2, 1.5)) {
  #   stopifnot( all.equal(1 - fGarch::.psnorm(q, xi), psnorm_upper_tail__std(q, xi), tolerance = 10^-12) )
  # }
}



#' MLE estimate of parameter Xi in a skew-normal distribution
#'
#' @param x numeric vector
#' @param mu mean value for x
#' @param sd sd value for x
#' @param init_xi initial value for xi (typically 1~1.25)
#' @param lower lower bound for Xi
#' @param upper Upper bound for Xi
fit_snorm_xi = function(x, mu, sd, init_xi, lower = NA, upper = NA) {

  ###### for reference; pure R variant (much slower than our Rcpp variant)
  # dsnorm_loglik = function(par) {
  #   -sum( dsnorm_log(x = x_standardized, xi = par) - log(sd) )
  # }
  # dsnorm_log = function(x, xi) {
  #   m1 = 2/sqrt(2 * pi)
  #   mu = m1 * (xi - 1/xi)
  #   g = 2/(xi + 1/xi)
  #   sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
  #   z = x * sigma + mu
  #   Xi = xi^sign(z)
  #   log(g) + dnorm(x = z/Xi, log = TRUE) + log(sigma)
  # }
  # x_standardized = (x - mu) / sd
  ######

  dsnorm_loglik = function(par) {
    rcpp_dsnorm_logsum(x_min_mean = x_min_mu, N = N, xi = par, sd = sd)
  }

  if(is.na(lower)) {
    lower = 1
  }
  if(is.na(upper)) {
    upper = 1.5
  }

  N = length(x)
  x_min_mu = x - mu

  stats::optim(
    par = init_xi, # initial parameters
    fn = dsnorm_loglik,
    method = "Brent",
    lower = lower,
    upper = upper,
    control = list(reltol = 1e-10)
  )
}
