
#' Wrapper function for goat_nulldistribution_function, typically used to generate the Rdata file with precomputed null distribution that is bundled with GOAT
#'
#' @examples \donttest{
#' # note that these example take a long time to compute
#'
#' # example 1; precompute parameters for a huge range of genesets, 500k permuations
#' goat_nulldistribution_function_generator(
#'   seq(from=100, to=20000, by=1), 500000, "C:/temp/", overwrite=TRUE, verbose=TRUE
#' )
#' # example 2; 2million permutations, a few very different genelist lengths; review QC plots
#' goat_nulldistribution_function_generator(
#'   c(100, 500, 1000, 5000, 15000), 2000000, "C:/temp/", overwrite=TRUE, verbose=TRUE
#' )
#' }
#' @param all_geneset_sizes vector of geneset sizes; integer values larger than 10 and smaller than half the genelist length (gene_scores length)
#' @param niter see `goat_nulldistribution_function`
#' @param output_dir full path to the directory where the downloaded files should be stored. Directory is created if it does not exist.
#' e.g. `output_dir="~/data"` on unix systems, `output_dir="C:/data"` on Windows, or set to `output_dir=getwd()` to write output to the current working directory
#' @param overwrite should existing files be overwritten? When setting this to FALSE, be mindful that you should manually clear cached files in between analyses !
#' @param verbose see `goat_nulldistribution_function`
#' @noRd
goat_nulldistribution_function_generator = function(all_geneset_sizes, niter, output_dir, overwrite = TRUE, verbose = FALSE) {

  geneset_size_bins = suppressWarnings(split(all_geneset_sizes, 1:min(100, length(all_geneset_sizes))))

  for(index_bin in seq_along(geneset_size_bins)) {
    start_time_bin = Sys.time()
    file_bin = sprintf("%s/goat_null_functions_%d.rda", output_dir, index_bin)
    if(overwrite == FALSE && file.exists(file_bin)) {
      next
    }

    result_genelists = list()
    for(genelist_N in geneset_size_bins[[index_bin]]) {
      start_time = Sys.time()
      genelist_fit = goat_nulldistribution_function(genelist_N, niter = niter, return_fit_objects = FALSE, verbose = (verbose && (genelist_N %in% c(100,250,500,750,1000,1500,2500,5000,10000,15000,20000))) )
      result_genelists[[length(result_genelists) + 1]] = genelist_fit
      if(index_bin == 1) {
        message(sprintf("N:%s took %d seconds", genelist_N, ceiling(as.numeric(difftime(Sys.time(), start_time, units = "secs"))) ))
      }
    }

    # store current batch to disk
    save(result_genelists, file = file_bin)
    message(sprintf("bin:%d took %.1f hours", index_bin, as.numeric(difftime(Sys.time(), start_time_bin, units = "hours")) ))
  }


  # from cached Rdata files per 'bin' into 1 large data.frame
  alldata = list()
  for(f in dir(output_dir, pattern = "^goat_null_functions_\\d+.rda$", full.names = TRUE)) {
    load(f)
    alldata[[length(alldata) + 1]] = result_genelists
  }
  goat_nulldistributions = dplyr::bind_rows(alldata) |> mutate_all(unname)

  # store final result
  save(goat_nulldistributions, file = paste0(output_dir, "/goat_null_alldata.rda"), compress = "xz", compression_level = 9)
}



#' Estimate null distribution parameters for a given genelist length as a polynomial function of geneset sizes
#'
#' 1) For a range of geneset sizes, from smallest (10) to largest (genelist length / 2), generate permutation-based null distributions and estimate their skew-normal parameters mean, sd, xi
#' 2) Fit monotonic polynomial functions across geneset sizes that predict sd or xi
#'
#' @param genelist_N genelist length. Integer value between 100 and 50000
#' @param niter number of bootstrap iterations that should be performed for generating empirical null distributions. Note that setting a large number will dramatically increase the RAM usage of this function! Integer value between 10000 and 5000000
#' @param return_fit_objects should the 'fit' object, which can be used as input for `predict`, be included in results? boolean value, TRUE or FALSE, default; FALSE
#' @param verbose should plots be created that describe the accuracy of fitting sd and xi? boolean value, TRUE or FALSE, default; FALSE
#' @noRd
goat_nulldistribution_function = function(genelist_N, niter = 500000, return_fit_objects = FALSE, verbose = FALSE) {
  validate_goat_niter(niter)
  validate_goat_verbose(verbose)
  stopifnot("parameter genelist_N (genelist length) should be an integer value >= 100 and <= 50000" =
              length(genelist_N) == 1 && is.numeric(genelist_N) && is.finite(genelist_N) && genelist_N >= 100 && genelist_N <= 50000)
  stopifnot("return_fit_objects parameter should be TRUE or FALSE" = length(return_fit_objects) == 1 && return_fit_objects %in% c(TRUE, FALSE))

  genelist_N = as.integer(genelist_N)
  niter = as.integer(niter)


  # compute an array of gene scores for the given genelist length
  # the order does not matter; this is the overall distribution of gene scores to sample from when computing geneset scores
  gene_scores = rankscore_fixed_order(genelist_N)


  #### DEFINE GENESET SIZES/BINS TO FIT
  max_geneset_size = min(10000, ceiling(genelist_N/2)) # hard cap to at most 10k genesets

  # approximately 250 bins with most 'density' around areas where we expect skew (and changes in skew between bins)
  if(max_geneset_size <= 255) {
    usizes = seq.int(from=5, to=max_geneset_size, by=1)
  } else if(max_geneset_size <= 500) {
    usizes = c(seq.int(from=5, to=180, by=1),
               seq.int(from=182, to=250, by=2),
               seq.int(from=255, to=max_geneset_size, by=5))
  } else if(max_geneset_size <= 1500) {
    usizes = c(seq.int(from=5, to=150, by=1),
               seq.int(from=152, to=200, by=2),
               seq.int(from=205, to=350, by=5),
               seq.int(from=360, to=490, by=10),
               seq.int(from=500, to=max_geneset_size, by=25))
  } else if(max_geneset_size <= 3000) {
    usizes = c(seq.int(from=5, to=150, by=1),
               seq.int(from=152, to=200, by=2),
               seq.int(from=205, to=250, by=5),
               seq.int(from=260, to=400, by=10),
               seq.int(from=245, to=750, by=25),
               seq.int(from=750, to=max_geneset_size, by=50))
  } else {
    usizes = c(seq.int(from=5, to=150, by=1),
               seq.int(from=152, to=200, by=2),
               seq.int(from=205, to=250, by=5),
               seq.int(from=260, to=350, by=10),
               seq.int(from=375, to=500, by=25),
               seq.int(from=550, to=1500, by=50),
               seq.int(from=1600, to=2500, by=100),
               seq.int(from=2750, to=max_geneset_size, by=250))
  }
  usizes = unique(c(usizes, max_geneset_size))


  #### SAMPLE NULL DISTRIBUTIONS
  mu = mean(gene_scores)
  null_data = rcpp_null_distributions(gene_scores = gene_scores, geneset_sizes = usizes, max_geneset_size = max_geneset_size, niter = niter)
  # init table with all stats per geneset bin/size
  df = data.frame(genelist_N = genelist_N, size = usizes, mu = mu, sigma = null_data$sigma, xi = 1)


  #### FIT SIGMA FUNCTION
  # fit a polynomial function to the observed standard deviations
  # Note that we include geneset sizes/bins larger than 50% of N, which won't be included in the results but are still useful to ensure we get a proper fit to the tail-end of the curve
  x = log(df$size)
  y = df$sigma
  flag = df$size >= 10 & df$size <= ceiling(genelist_N/2)
  x = x[flag]
  y = y[flag]
  # fit using MonoPoly package
  k_sigma = 9
  fit_sigma = MonoPoly::monpol(y ~ x, degree = k_sigma)
  df$sigma_smooth = as.numeric(stats::predict(fit_sigma, newdata = data.frame(x = log(df$size))))
  if(verbose) {
    # compute residual standard error, ignoring weights. To include weights, use; summary(fit)$sigma
    fit_se = sqrt(sum(fit_sigma$residuals^2) / (length(fit_sigma$residuals) - length(fit_sigma$coefficients)))
    plot(x, y, cex = 0.5, pch = c(16,1)[1 + (exp(x) < 10 | exp(x) > genelist_N/2)], xlab = "log(size)", ylab = "sd", main = sprintf("genelist N:%d  Sigma residual standard error: %.5g", genelist_N, fit_se), cex.main = 0.9)
    graphics::lines(x, stats::predict(fit_sigma), col = 2)
  }


  #### FIT XI AT RAW NULL DATA PER BIN/SIZE
  prev_xi = NA # set NA to indicate that we haven't fitted any actual data yet
  for(index in 1L:nrow(df)) {
    x = null_data$null[(1 + (index-1) * niter):(index * niter)]
    # note the use of smoothed SD estimates + limit Xi possible params to "near previous estimate" for a minor speedup of the MLE fit
    iter_mle_fit = fit_snorm_xi(x, mu, sd = df$sigma_smooth[index], init_xi = ifelse(is.na(prev_xi), 1.25, prev_xi), lower = ifelse(is.na(prev_xi), 1, max(1, prev_xi * 0.975)), upper = ifelse(is.na(prev_xi), 1.5, prev_xi * 1.025) )
    df$xi[index] = prev_xi = iter_mle_fit$par
  }


  #### FIT XI FUNCTION
  x = log(df$size)
  y = log(df$xi)
  flag = df$size >= 10 & df$size <= ceiling(genelist_N/2)
  x = x[flag]
  y = y[flag]
  # fit using MonoPoly package
  k_xi = 11
  fit_xi = MonoPoly::monpol(y ~ x, degree = k_xi, weights = df$xi[flag] - 0.5) # weight = increase importance of increasingly-skewed geneset bins
  # df$xi_smooth = exp(as.numeric(stats::predict(fit_xi, data.frame(x=log(df$size))))) # only needed for debugging
  if(verbose) {
    # compute residual standard error, ignoring weights. To include weights, use; summary(fit)$sigma
    fit_se = sqrt(sum(fit_xi$residuals^2) / (length(fit_xi$residuals) - length(fit_xi$coefficients)))
    plot(x, y, cex = 0.5, pch = c(16,1)[1 + (exp(x) < 10 | exp(x) > genelist_N/2)], xlab = "log(size)", ylab = "log(Xi)", main = sprintf("genelist N:%d  Xi residual standard error: %.5g", genelist_N, fit_se), cex.main = 0.9)
    graphics::lines(x, stats::predict(fit_xi), col = 2)
  }


  #### store final result
  ## QC output from fit functions just to be sure
  coef_sigma = stats::coef(fit_sigma) # intercept followed by coefficients of x^k
  coef_xi    = stats::coef(fit_xi)
  stopifnot(length(coef_sigma) == k_sigma + 1 && all(is.finite(coef_sigma)))
  stopifnot(length(coef_xi)    == k_xi + 1    && all(is.finite(coef_xi)))
  coef_sigma = as.numeric(coef_sigma) # enforce plain numeric vector as type
  coef_xi    = as.numeric(coef_xi)

  ## store data in list format with named parameters for Sigma and Xi  (where k=0 is the intercept and others are coef*x^k)
  genelist_fit = list(N=genelist_N, mu=mu)
  for(i in seq_along(coef_sigma)) {
    genelist_fit[[paste0("sigma_", i-1)]] = coef_sigma[i]
  }
  for(i in seq_along(coef_xi)) {
    genelist_fit[[paste0("xi_", i-1)]] = coef_xi[i]
  }

  if(return_fit_objects) {
    genelist_fit$fit_sigma = fit_sigma
    genelist_fit$fit_xi = fit_xi
  }

  return(genelist_fit)
}



#' Estimate null distribution parameters for a given set of geneset sizes and an array of gene scores (from genelist)
#'
#' @param geneset_sizes vector of geneset sizes; integer values larger than 10 and smaller than half the genelist length (gene_scores length)
#' @param gene_scores gene score vector; these values are assumed to be the 'universe' of all possible gene-level scores. These are used to compute geneset score null distributions
#' @param niter number of bootstrap iterations that should be performed for generating empirical null distributions. Note that setting a large number will dramatically increase the RAM usage of this function! Integer value between 10000 and 5000000
#' @param verbose boolean, plot 10 null distribution histograms augmented with fitted skew-normal distributions
#' @return data.frame with columns; genelist_N, size, mu, sigma, xi
#' @noRd
goat_nulldistribution_independent = function(geneset_sizes, gene_scores, niter, verbose) {
  validate_goat_niter(niter)
  validate_goat_verbose(verbose)
  stopifnot(is.numeric(geneset_sizes) & is.finite(geneset_sizes))
  stopifnot(is.numeric(gene_scores) & is.finite(gene_scores) & gene_scores >= 0)

  niter = as.integer(niter)
  genelist_N = length(gene_scores)
  usizes = sort(unique(geneset_sizes))
  min_geneset_size = usizes[1]
  max_geneset_size = utils::tail(usizes, 1)
  stopifnot(min_geneset_size >= 10 && max_geneset_size <= min(10000, ceiling(genelist_N/2)))
  mu = mean(gene_scores)

  ### generate null distributions and compute sd() for each unique geneset size (in input data)
  # this C++ function will use multiprocessing if openmp support is available
  null_data = rcpp_null_distributions(gene_scores = gene_scores, geneset_sizes = usizes, max_geneset_size = max_geneset_size, niter = niter)


  ### iterate geneset size bins and fit skew normal distributions. After converging to Gaussian at geneset size N, we stop fitting skew normal
  df = data.frame(genelist_N = genelist_N, size = usizes, mu = mu, sigma = null_data$sigma, xi = 1)
  prev_xi = 1.15
  for(index in 1L:nrow(df)) {
    # from the null distribution data, an array that encodes the "geneset_size * bootstrap_iteration" matrix by row,
    # we extract the null distribution for current geneset size bin = usizes[index]
    x = null_data$null[(1 + (index-1) * niter):(index * niter)]

    iter_mle_fit = fit_snorm_xi(x, mu, sd = df$sigma[index], init_xi = prev_xi, lower = 1, upper = 1.5)
    # if skew-normal, store result and continue to next geneset size bin
    if(iter_mle_fit$par >= 1.005) {
      df$xi[index] = prev_xi = iter_mle_fit$par
    } else {
      break
    }
  }


  ### QC; plot raw data versus fits
  if(verbose) {
    ## find at most 10 geneset sizes that are interesting for QC plot; a few of the smallest + bins closest to 25, 50, 100, 150, 250, 500
    if(nrow(df) <= 10) {
      index_plot = seq_len(nrow(df))
    } else {
      # add smallest 2
      index_plot = 1L:2L
      # add at most 2 that are larger, but only by 3~5
      for(i in seq_len(2)) {
        j = which( df$size >= 3L + df$size[utils::tail(index_plot,1)] & df$size <= 5L + df$size[utils::tail(index_plot,1)] )
        if(length(j) > 0) {
          index_plot = c(index_plot, j[1])
        }
      }
      # now add closest to value-of-interest
      for(bin in c(25, 50, 100, 150, 250, 500)) {
        j = which(abs(df$size - bin) == min(abs(df$size - bin)))
        if(length(j) > 0) {
          index_plot = c(index_plot, j[1])
        }
      }

      # done. enforce uniqueness
      index_plot = unique(index_plot)
    }


    # iterate bins (index in df) and plot
    for(index in index_plot) {
      # analogous to above df loop, grab null distribution data for geneset of current size
      x = null_data$null[(1 + (index-1) * niter):(index * niter)]
      x_qtl = stats::quantile(x, probs = c(0.0001, 0.9999))
      x = x[x > x_qtl[1] & x < x_qtl[2]] # limit data to plot range
      # plot histogram
      graphics::hist(x, breaks = 90, freq = FALSE, xlab = "geneset score", ylab = "density", main = sprintf(
        "null distribution for geneset size:%d = sd:%.3f xi:%.3f", df$size[index], df$sigma[index], df$xi[index]
      ), cex.main = 0.9)
      # add a line for the fitted null distribution (from which we estimate the actual geneset p-values later on)
      lgnd = c("normal", "skew-normal fit")
      lgnd_col = c("blue", "red")
      graphics::curve(stats::dnorm(x, mean = mu, sd = df$sigma[index]), add = TRUE, col = "blue", lwd = 2)
      graphics::curve(dsnorm(x, mean = mu, sd = df$sigma[index], xi = df$xi[index]), add = TRUE, col = "red", lwd = 2, lty = 2)
      graphics::legend("topright", legend = lgnd, col = lgnd_col, lwd = 2)
    }

  }

  return(df)
}
