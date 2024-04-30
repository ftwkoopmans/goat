

# this function is not integrated into the GOAT R package atm as it causes issues while submitting the GOAT R package to CRAN


#' Geneset testing with the iDEA R package
#'
#' @description
#' Follows the iDEA tutorial at https://xzhoulab.github.io/iDEA/ with minor adaptions (see `idea_variant` parameter).
#'
#' note that computation takes a long time with the iDEA method. On a high performance workstation computer with
#' 12 cores / 24 threads it takes approximately 6 hours for a genelist of ~14000 genes to analyze all (filtered) GO
#' genesets. Using a regular PC this will take much longer.
#'
#'
#' @param genesets filtered genesets
#' @param genelist user genelist, must contain columns "gene", "log2fc" and "pvalue"
#' @param return_pvalue_louis boolean, set `TRUE` (default) to return "louis corrected" p-values from iDEA (analogous to iDEA tutorial)
#' @param idea_variant single string indicating how iDEA should be used. Valid options: 'default', 'beta', 'beta_rescale', 'beta_rescale_center', 'rescale', 'rescale_center', 'genelist_sd', 'genelist_var'.
#' default = follow iDEA tutorial
#' beta = let iDEA model only the beta's (log2fc), this is equivalent to `modelVariant = TRUE` in the iDEA tutorial
#' beta_rescale = analogous to 'beta', but with beta's rescaled by their MAD (i.e. we ensure that iDEA receives input data at expected scale)
#' beta_rescale_center = analogous to 'beta_rescale', but also center to mode zero
#' rescale = rescale input such that beta_var is approximately gamma(2,0.5) distributed. This resolves the vast majority of iDEA model fitting errors that otherwise occur at default settings
#' rescale_center = analogous, but the rescaled beta/log2fc are also centered to mode zero
#' genelist_sd = use provided standard deviations (of the log2fc) in the 'sd' column of the genelist
#' genelist_var = use provided variations (of the log2fc) in the 'var' column of the genelist
#' @param verbose boolean, set `TRUE` to plot genelist beta and beta/sqrt(beta_var) distributions
#' @returns input `genesets` with additional columns for geneset p-value and score type. Because iDEA can return 2 types of geneset p-values,
#' multiple columns are included in this function's results; `pvalue_idea` = iDEA geneset p-value, `pvalue_louis_idea` = `pvalue_louis` from iDEA,
#' `pvalue` = respective copy of either former column, depending on user-parameter `return_pvalue_louis`
#' @seealso `test_genesets`
#' @return input `genesets` table with results in the "pvalue" and "score_type" columns.
#' Additional data at "ismissing_idea", "pvalue_idea" and "pvalue_louis_idea" columns
#' @export
test_genesets_idea = function(genesets, genelist, return_pvalue_louis = TRUE, idea_variant = "default", verbose = FALSE) {
  x = NULL # fix invisible bindings R package NOTE
  force(genesets) # guard against issues with multithreaded code downstream; enforce evaluation of input data (i.e. no lazy evaluation)
  force(genelist)
  if(!requireNamespace("future.apply", quietly = TRUE)) {
    stop('the "future.apply" R package is required for test_genesets_idea(). Start a fresh RStudio session and issue this R statement: install.packages("future.apply)')
  }
  if(!requireNamespace("future", quietly = TRUE)) {
    stop('the "future" R package is required for test_genesets_idea(). Start a fresh RStudio session and issue this R statement: install.packages("future)')
  }
  if(!requireNamespace("iDEA", quietly = TRUE)) {
    stop('to perform geneset testing with the iDEA algorithm, the iDEA R package is required which is available on GitHub.\nSuggested installation steps;\n1) close RStudio and start it anew (do copy/paste these instructions before closing)\n2) run R command: install.packages("pak")\n3) run R command: pak::pkg_install("xzhoulab/iDEA", update=FALSE)\n4) test if the installation was successful by loading the R package, run R command: library(iDEA)', call. = FALSE)
  }


  # input validation
  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns; source (optional, character type) id (character type), genes (list type)" =
              is.data.frame(genesets) && nrow(genesets) > 0 && all(c("id", "genes") %in% colnames(genesets)) &&
              is.character(genesets$id) && is.list(genesets$genes))
  stopifnot("genelist parameter must be a non-empty data.frame/tibble with numeric columns that contain finite values; 'gene' (integer), 'log2fc' (numeric) and 'pvalue' (numeric)" =
              is.data.frame(genelist) && nrow(genelist) > 0 && all(c("gene", "log2fc", "pvalue") %in% colnames(genelist)) &&
              all(is.finite(genelist$gene)) && all(is.finite(genelist$log2fc)) && all(is.finite(genelist$pvalue)))
  stopifnot("return_pvalue_louis parameter must be TRUE or FALSE" = length(return_pvalue_louis) == 1 && return_pvalue_louis %in% c(TRUE, FALSE))
  stopifnot("verbose parameter must be TRUE or FALSE" = length(verbose) == 1 && verbose %in% c(TRUE, FALSE))
  stopifnot("idea_variant parameter should be a single string. Valid options: 'default', 'beta', 'beta_rescale', 'beta_rescale_center', 'rescale', 'rescale_center', 'genelist_sd', 'genelist_var'" = length(idea_variant) == 1 && idea_variant %in% c("default", "beta", "beta_rescale", "beta_rescale_center", "rescale", "rescale_center", "genelist_sd", "genelist_var"))
  if(idea_variant == "genelist_sd") {
    stopifnot("idea_variant parameter is genelist_sd, expect a column named 'sd' in provided genelist with only finite numeric values larger than 0" = "sd" %in% colnames(genelist) && all(is.numeric(genelist$sd) & is.finite(genelist$sd) & genelist$sd > 0))
  }
  if(idea_variant == "genelist_var") {
    stopifnot("idea_variant parameter is genelist_var, expect a column named 'var' in provided genelist with only finite numeric values larger than 0" = "var" %in% colnames(genelist) && all(is.numeric(genelist$var) & is.finite(genelist$var) & genelist$var > 0))
  }

  # init results
  genesets$score_type = NA_character_
  genesets$pvalue = NA_real_

  # temp geneset ID that combines source + id for compatibility with input genelists that contain the same ID across 'source'
  genesets$id_tmp = genesets$id
  if("source" %in% colnames(genesets)) {
    genesets$id_tmp = paste(genesets$source, genesets$id)
  }


  ### genelist to iDEA format (data.frame with iDEA-required 'beta' and 'beta_var' columns)
  genelist_idea = idea_prepare_genelist(genelist, idea_variant, verbose)
  # this is equivalent to iDEA parameter "modelVariant".
  # From iDEA tutorial; "model option to run, boolean variable, if FALSE, runing the main iDEA model, which models on z score statistics. if TRUE, runing iDEA variant model which models on beta effect size."
  idea_param_ismodelvariant = idea_variant %in% c("beta", "beta_rescale", "beta_rescale_center")
  cat("iDEA variant: ", idea_variant, "\n", sep = "")


  ### genesets to iDEA format (gene*geneset identity matrix)
  # helper function creates the identity matrix, which we can later trivially convert to a data.frame if needed
  genesets_to_matrix = function(gs, ugene) {
    # note; in GOAT we use numeric gene IDs. To avoid type conversion / indexing shenanigans, we prefix integer gene IDs with a string
    m = matrix(0L, nrow = length(ugene), ncol = nrow(gs), dimnames = list(paste0("g", ugene), gs$id_tmp)) # ! use the temp ID for genesets
    for(i in 1:nrow(gs)) {
      m[,i] = as.integer(ugene %in% gs$genes[[i]])
    }
    return(m)
  }
  genesets_idea = as.data.frame(genesets_to_matrix(genesets, genelist$gene))
  # double-check data structure integrity
  stopifnot(rownames(genesets_idea) == rownames(genelist_idea))


  ### following the iDEA tutorial, perform geneset enrichment testing with default parameters
  # setup multiprocessing; start future plan + cleanup on exit
  oplan = future::plan(future::multisession)
  on.exit(future::plan(oplan))

  tryCatch({
    # as first parameter, select only the 2 required input columns and enforce their order (i.e. this is a hardcoded requirement in iDEA)
    # we set `min_precent_annot=0` because we already filtered genesets upstream (default in tutorial is 0.0025, which amounts to 35 in a genelist of 14k genes; skipping over small genesets)
    # we disabled multiprocessing here because it has issues on some of our test systems (Windows); this step is quite fast so no problem
    idea = iDEA::CreateiDEAObject(genelist_idea[,c("beta", "beta_var")], genesets_idea, max_var_beta = 101, min_precent_annot = 0, num_core = 1)

    ################
    # idea = iDEA::iDEA.fit(
    #   idea,
    #   fit_noGS = FALSE,
    #   init_beta = NULL,
    #   init_tau = c(-2,0.5),
    #   min_degene = 0,
    #   em_iter = 15,
    #   mcmc_iter = 1000,
    #   fit.tol = 1e-5,
    #   modelVariant = F,
    #   verbose = TRUE
    # )
    ################

    ### ANALOGOUS TO ABOVE CODE FROM IDEA TUTORIAL, DIRECTLY ADAPTED FROM iDEA SOURCE CODE, BUT WITH future.apply TO FIX PARALLEL COMPUTATION  +  SWITCH VARIANT @ EMMCMCStepSummaryVariant / EMMCMCStepSummary
    # https://github.com/xzhoulab/iDEA/blob/a30fc8ec655813678049c5b4b4b0aa6a21c944f0/R/iDEASummary.R#L103
    # here setting default parameters from tutorial at: https://xzhoulab.github.io/iDEA/
    object = idea
    init_tau = c(-2,0.5) # initial value for the coefficient of annotations/gene sets, including the intercept in EM procedure, default is c(-2,0.5).
    em_iter = 15         # maximum iteration for EM algorithm, default is 15
    mcmc_iter = 1000     # maximum iteration for MCMC algorithm, default is 1000
    min_degene = 0       # the threshold for the number of detected DE genes in summary statistics. For some of extremely cases, the method does not work stably when the number of detected DE genes is 0.
    #### start snippet from iDEA source
    num_gene <- object@num_gene
    num_annot <- length(object@annotation)
    init_beta <- object@summary[,1]
    # this line is the only deviation from the iDEA sourcecode;
    # we use the future.apply package because the original code, using 'pbmclapply'
    # 1) occasionally hangs in WSL and 2) doesn't perform well on Windows
    # (this is a simple workaround as we only use this specific iDEA configuration)
    # (note that we can't print to console from within future "workers", so iDEA-specific-errors are lost / not reported)
    res_idea <- future.apply::future_sapply(1:num_annot, FUN = function(x) {
      set.seed(x)
      Annot <- rep(0, object@num_gene)
      Annot[object@annotation[[x]]] <- 1
      Annot <- Annot - mean(Annot)
      Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
      if(idea_param_ismodelvariant) {
        t1 <- system.time(model1 <- try( res <- iDEA::EMMCMCStepSummaryVariant(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
      } else {
        t1 <- system.time(model1 <- try( res <- iDEA::EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
      }
      if(!inherits(model1, "try-error")) {
        rownames(res$pip) <- object@gene_id
        colnames(res$pip) <- "PIP"
        rownames(res$beta) <- object@gene_id
        colnames(res$beta) <- "beta"
        rownames(res$annot_coef) <- c("tau_1", "tau_2")
        colnames(res$annot_coef) <- "annot_coef"
        rownames(res$annot_var) <- c("tau_1", "tau_2")
        colnames(res$annot_var) <- "annot_var"
        res$converged   <- TRUE
        res$ctime   <- t1[3]
      } else { res <- NULL }
      return(res)
    }, future.seed = NULL, simplify = FALSE)

    names(res_idea) <- object@annot_id
    object@de <- res_idea
    #### end snippet from iDEA source
    idea@de <- res_idea
    ################

    ### fix the the estimated geneset p-values as per default iDEA tutorial workflow
    tmp = utils::capture.output(idea <- iDEA::iDEA.louis(idea)) # capture.output to silence the progress bar

    ### map results to input geneset data.frame, from the idea object
    index = match(genesets$id_tmp, idea@gsea$annot_id) # ! use the temp ID for genesets
    if(any(!is.finite(idea@gsea$pvalue_louis[index]))) {
      cat(sprintf("%d/%d (filtered) genesets that were tested with iDEA received a p-value via iDEA.louis(), iDEA model fitting failed for %d genesets\n",
                  sum(is.finite(idea@gsea$pvalue_louis[index])), nrow(genesets), sum(!is.finite(idea@gsea$pvalue_louis[index])) ))
    }
    genesets$ismissing_idea = !is.finite(idea@gsea$pvalue_louis[index])
    genesets$pvalue_idea = idea@gsea$pvalue[index]
    genesets$pvalue_louis_idea = idea@gsea$pvalue_louis[index]
    genesets$pvalue = idea@gsea$pvalue[index]
    if(return_pvalue_louis) {
      genesets$pvalue = idea@gsea$pvalue_louis[index]
    }
    # In this function we tested only valid, non-empty genesets that passed filter_genesets() upstream.
    # If the geneset algorithm returned a non-finite pvalue, we consider it "tested but found not significant"
    # (note that leaving NA values as-is would affect subsequent multiple testing correction, which in turn would be a
    # source of bias if many "no effect whatsoever" genesets were returned as NA pvalues by the algorithm)
    genesets$pvalue_idea[!is.finite(genesets$pvalue_idea)] = 1
    genesets$pvalue_louis_idea[!is.finite(genesets$pvalue_louis_idea)] = 1
    genesets$pvalue[!is.finite(genesets$pvalue)] = 1
    # score type
    genesets$score_type = c("effectsize_down", "effectsize_up")[1 + (is.finite(idea@gsea$annot_coef[index]) & idea@gsea$annot_coef[index] > 0)] # if annot_coef > 0, pathway is up-regulated
  },
  error = function(err) {
    print(err)
  })

  # remove temp geneset ID column
  genesets$id_tmp = NULL
  return(genesets)
}



#' prepare genelist table with 'beta' and 'beta_var' variables depending on desired iDEA variant (e.g. only use beta and/or rescale input values)
#'
#' @description
#' ##### for reference, relevant snippet from iDEA tutorial @ https://xzhoulab.github.io/iDEA/ that we adopted here
#' `## Assume you have obtained the DE results from i.e. zingeR, edgeR or MAST with the data frame res_DE (column: pvalue and LogFC)`
#' `pvalue <- res_DE$pvalue #### the pvalue column`
#' `zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score`
#' `beta <- res_DE$LogFC ## effect size`
#' `se_beta <- abs(beta/zscore) ## to approximate the standard error of beta`
#' `beta_var = se_beta^2  ### square`
#' `summary = data.frame(beta = beta,beta_var = beta_var)`
#' `## add the gene names as the rownames of summary`
#' `rownames(summary) = rownames(res_DE) ### or the gene id column in the res_DE results`
#'
#' @param genelist see `test_genesets_idea`
#' @param idea_variant see `test_genesets_idea`
#' @param verbose see `test_genesets_idea`
#' @noRd
idea_prepare_genelist = function(genelist, idea_variant, verbose) {
  genelist_idea = x = NULL

  ### genelist to iDEA format (assuming we have columns log2fc and pvalue @ input genelist)
  if(idea_variant == "genelist_sd") {
    ### user-provided input values
    genelist_idea = data.frame(beta = genelist$log2fc, beta_var = genelist$sd^2, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
  }
  if(idea_variant == "genelist_var") {
    ### user-provided input values
    genelist_idea = data.frame(beta = genelist$log2fc, beta_var = genelist$var, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
  }

  # rescale beta up-front @ genelist, then follow typical iDEA genelist preparation
  if(idea_variant %in% c("beta_rescale", "beta_rescale_center")) {
    # rescaling factor is basically the log2fc MAD. Below implementation increases robustness by ignoring zeros and extreme values
    b_threshold = 3 * stats::quantile(abs(genelist$log2fc[is.finite(genelist$log2fc) & genelist$log2fc != 0]), probs = 0.9)
    k = stats::mad(genelist$log2fc[is.finite(genelist$log2fc) & genelist$log2fc != 0 & abs(genelist$log2fc) < b_threshold])
    genelist$log2fc = genelist$log2fc / k
  }
  # note; in GOAT we use numeric gene IDs. To avoid type conversion / indexing shenanigans, we prefix integer gene IDs with a string
  if(idea_variant %in% c("default", "beta", "beta_rescale", "beta_rescale_center")) {
    ### here we basically follow the iDEA tutorial, but with minor adaptions for edge-cases that cause zero or infinite se_beta estimates (e.g. log2fc=0 or pvalue=1)
    genelist_idea = data.frame(beta = genelist$log2fc, beta_var = NA, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
    genelist_idea$zscore = stats::qnorm(genelist_idea$pvalue / 2.0, lower.tail=FALSE)
    genelist_idea$se_beta = abs(genelist_idea$beta / genelist_idea$zscore) # approximate the standard error of beta (log2fc)
    genelist_idea$se_ismissing = !is.finite(genelist_idea$se_beta) | genelist_idea$se_beta == 0
    genelist_idea$se_beta[genelist_idea$se_ismissing] = max(genelist_idea$se_beta[!genelist_idea$se_ismissing]) # update/fix; insert high values where beta_se is missing
    genelist_idea$se_beta[genelist_idea$se_beta < 10^-6] = 10^-6 # update/fix; guard against "zero standard error" (e.g. at log2fc=0; 0/x=0)
    genelist_idea$beta_var = genelist_idea$se_beta^2
  }

  if(idea_variant %in% c("rescale", "rescale_center")) {
    ### follow the iDEA tutorial, but additionally we standardize beta and beta_var
    ## first, we follow the iDEA tutorial @ https://xzhoulab.github.io/iDEA/
    genelist_idea = data.frame(beta = genelist$log2fc, beta_var = NA, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
    genelist_idea$zscore = stats::qnorm(genelist_idea$pvalue / 2.0, lower.tail=FALSE)
    genelist_idea$se_beta = abs(genelist_idea$beta / genelist_idea$zscore) # approximate the standard error of beta (log2fc)
    genelist_idea$se_ismissing = !is.finite(genelist_idea$se_beta) | genelist_idea$se_beta == 0
    genelist_idea$se_beta[genelist_idea$se_ismissing] = max(genelist_idea$se_beta[!genelist_idea$se_ismissing]) # update/fix; insert high values where beta_se is missing
    genelist_idea$beta_var = genelist_idea$se_beta^2

    ## next; rescale such that beta_var is always on a similar scale regardless of the input data
    # This will fix the many iDEA errors observed in datasets that have atypical SD distributions (when using the 'default' iDEA workflow).
    # The simple scaling factor will result in roughly gamma(2, 0.5) shaped distributions
    # TODO:
    # we could make this more robust with a basic MLE to rescale beta_var distribution to expectation from gamma(2, 0.5),
    # instead of rescaling to a constant as used here (but for a dozen datasets we tested, this works reasonably well and at least fixes all iDEA errors)
    # Analogous to the threshold/non-zero filtering prior to MAD for beta_rescale, this also ignores extreme outliers (and bimodality when there are many missing SE that were capped at extreme/max value)
    bv_threshold = 2 * stats::quantile(genelist_idea$beta_var[!genelist_idea$se_ismissing], probs = 0.9)
    k = 0.6 / stats::mad(genelist_idea$beta_var[!genelist_idea$se_ismissing & genelist_idea$beta_var < bv_threshold])
    genelist_idea$std__var = genelist_idea$beta_var * k
    genelist_idea$std__var[genelist_idea$std__var < 10^-6] = 10^-6 # guard against "zero standard error" (e.g. caused by log2fc=0)
    genelist_idea$std__var[genelist_idea$std__var > 50] = 50 # after standardizing, we can cap extreme values
    genelist_idea$std__beta = genelist_idea$beta * sqrt(k) # rescale beta as well, so beta/sd ratio remains unchanged
    # result: new table with only the updated values
    genelist_idea = data.frame(beta = genelist_idea$std__beta, beta_var = genelist_idea$std__var, pvalue = genelist_idea$pvalue, row.names = paste0("g", genelist$gene))
  }

  ## optionally, center the foldchanges to zero mode
  if(idea_variant %in% c("rescale_center", "beta_rescale_center")) {
    get_mode = function(x) {
      density_estimate = stats::density(x, na.rm=T)
      return(density_estimate$x[which.max(density_estimate$y)])
    }
    m = get_mode(genelist_idea$beta) # in contrast to above code for rescaling, here we use all log2fc values (including zero's and those with pvalue=1)
    stopifnot("failed to estimate mode of (centered) beta @ iDEA" = is.finite(m))
    genelist_idea$beta = genelist_idea$beta - m
  }

  if(is.null(genelist_idea)) {
    stop(paste("unknown option for iDEA;", idea_variant))
  }


  if(verbose) {
    mylim = function(x, qtl = c(0.0025, 0.9975)) c(-1, 1) * max(abs(stats::quantile(x[is.finite(x) & x != 0], probs = qtl, na.rm = TRUE)))
    mydens = function(x, lim, adj = 0.5) stats::density(x[is.finite(x) & x >= lim[1] & x <= lim[2]], adjust = adj)
    oldpar = graphics::par(mfcol = c(3,1))

    d = mydens(genelist_idea$beta, mylim(genelist_idea$beta))
    plot(d, xlab = "beta", main = paste0("iDEA variant: ", idea_variant, "\ngrey line: standard normal"))
    graphics::abline(v = 0, col = "grey")
    graphics::curve(stats::dnorm(x, mean = 0, sd = 1), add = TRUE, col = "grey")

    d = mydens(genelist_idea$beta_var, c(0, stats::quantile(genelist_idea$beta_var, probs = 0.9975)) )
    plot(d, xlab = "variation", main = "grey line: gamma(shape=2, scale=0.5)")
    graphics::curve(stats::dgamma(x, shape = 2, scale = 0.5), add = TRUE, col = "grey")

    d = mydens(genelist_idea$beta / sqrt(genelist_idea$beta_var), mylim(genelist_idea$beta / sqrt(genelist_idea$beta_var)))
    plot(d, xlab = "beta / sd", main = "grey line: standard normal")
    graphics::abline(v = 0, col = "grey")
    graphics::curve(stats::dnorm(x, mean = 0, sd = 1), add = TRUE, col = "grey")
    graphics::par(oldpar)
  }

  return(genelist_idea)
}
