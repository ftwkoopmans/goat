
#' Geneset testing with the iDEA R package (closely follows the iDEA tutorial)
#'
#' @description
#' note that computation takes a long time with the iDEA method. On a high performance workstation computer with
#' 12 cores / 24 threads it takes approximately 6 hours for a genelist of ~14000 genes to analyze all (filtered) GO
#' genesets. Using a regular PC this might take a multitude of hours.
#'
#' @param genesets filtered genesets
#' @param genelist user genelist, must contain columns "gene", "log2fc" and "pvalue"
#' @param return_pvalue_louis boolean, set `TRUE` to return "louis corrected" p-values from iDEA (analogous to iDEA tutorial)
#' @returns input `genesets` with additional columns for geneset p-value and score type. Because iDEA can return 2 types of geneset p-values,
#' multiple columns are included in this function's results; `pvalue_idea` = iDEA geneset p-value, `pvalue_louis_idea` = `pvalue_louis` from iDEA,
#' `pvalue` = respective copy of either former column, depending on user-parameter `return_pvalue_louis`
#' @seealso `test_genesets`
#' @export
test_genesets_idea = function(genesets, genelist, return_pvalue_louis) {
  force(genesets) # guard against issues with multithreaded code downstream; from lazy evaluation (default) to enforced
  force(genelist)
  check_dependency("iDEA", "geneset enrichment testing with iDEA")
  check_dependency("future.apply", "geneset enrichment testing with iDEA")
  # input validation
  stopifnot("genesets parameter must be a non-empty data.frame/tibble with columns; source (optional, character type) id (character type), genes (list type)" =
              is.data.frame(genesets) && nrow(genesets) > 0 && all(c("id", "genes") %in% colnames(genesets)) &&
              is.character(genesets$id) && is.list(genesets$genes))
  stopifnot("genelist parameter must be a non-empty data.frame/tibble with numeric columns that contain finite values; 'gene' (integer), 'log2fc' (numeric) and 'pvalue' (numeric)" =
              is.data.frame(genelist) && nrow(genelist) > 0 && all(c("gene", "log2fc", "pvalue") %in% colnames(genelist)) &&
              all(is.finite(genelist$gene)) && all(is.finite(genelist$log2fc)) && all(is.finite(genelist$pvalue)))
  stopifnot("return_pvalue_louis parameter must be TRUE or FALSE" = length(return_pvalue_louis) == 1 && return_pvalue_louis %in% c(TRUE, FALSE))

  # init results
  genesets$score_type = NA_character_
  genesets$pvalue = NA_real_

  # temp geneset ID that combines source + id for compatibility with input genelists that contain the same ID across 'source'
  genesets$id_tmp = genesets$id
  if("source" %in% colnames(genesets)) {
    genesets$id_tmp = paste(genesets$source, genesets$id)
  }


  ### genelist to iDEA format (assuming we have columns log2fc and pvalue @ input genelist)
  # this code generally follows the iDEA tutorial, but with minor adaptions for edge-cases that cause zero or infinite se_beta estimates
  # note; in GOAT we use numeric gene IDs. To avoid type conversion / indexing shenanigans, we prefix integer gene IDs with a string
  genelist_idea = data.frame(beta = genelist$log2fc, beta_var = NA, pvalue = genelist$pvalue, row.names = paste0("g", genelist$gene))
  # limit pvalues between 10^-64 and 1. also deals with upstream data that might (due to rounding errors or whatever in user input) input pvalue=0
  genelist_idea$zscore = stats::qnorm(pmin(1, pmax(genelist_idea$pvalue, 10^-64)) / 2.0, lower.tail = FALSE) # between ~17 and 0 (zcore=0 @ pvalue=1)
  genelist_idea$se_beta = abs(genelist_idea$beta / genelist_idea$zscore) # approximate the standard error of beta (log2fc)
  genelist_idea$se_beta[!is.finite(genelist_idea$se_beta)] = 10 # default to high value where beta_se is missing/invalid (zcore=0 @ pvalue=1 -> beta/0 = Inf or -Inf), yields beta_var = 100
  genelist_idea$se_beta[genelist_idea$se_beta < 10^-16] = 10^-16 # guard against "zero standard error" (e.g. at log2fc=0; 0/x=0) to ensure downstream iDEA code doesn't suffer from division-by-zero issues
  genelist_idea$beta_var = genelist_idea$se_beta^2


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
  ## alternatively, to completely remove genes that have no annotations from the genelist and genesets;
  # genesets_idea = as.data.frame(genesets_to_matrix(genesets, ugene = sort(unique(unlist(gs$genes, recursive = F, use.names = F)))))
  # genelist_idea = genelist_idea[match(rownames(genesets_idea), rownames(genelist_idea)), ]
  stopifnot(rownames(genesets_idea) == rownames(genelist_idea))

  # setup multiprocessing; start future plan + cleanup on exit
  oplan = future::plan(future::multisession)
  on.exit(future::plan(oplan))

  tryCatch({
    ### following the iDEA tutorial, perform geneset enrichment testing with default parameters
    # as first parameter, select only the 2 required input columns and enforce their order (i.e. this is a hardcoded requirement in iDEA)
    # we set `min_precent_annot=0` because we already filtered genesets upstream (i.e. default in tutorial is 0.0025, which amounts to 35 in a genelist of 14k genes; skipping over small genesets)
    # we disabled multiprocessing here because it has issues on some of our test systems; this step is quite fast so no problem
    idea = iDEA::CreateiDEAObject(genelist_idea[,c("beta", "beta_var")], genesets_idea, max_var_beta = 101, min_precent_annot = 0, num_core = 1)


    ################
    # idea = iDEA::iDEA.fit(
    #   idea,
    #   fit_noGS = FALSE,
    #   init_beta = NULL,
    #   init_tau = c(-2,0.5),
    #   min_degene = 0, # defaults to 5 in the iDEA tutorial; "min_degene: the threshold for the number of detected DE genes in summary statistics. For some of extremely cases, the method does not work stably when the number of detected DE genes is 0"
    #   em_iter = 15,
    #   mcmc_iter = 1000,
    #   fit.tol = 1e-5,
    #   modelVariant = F,
    #   verbose = TRUE
    # )
    ################  ANALOGOUS TO ABOVE CODE FROM IDEA TUTORIAL, DIRECTLY ADAPTED FROM iDEA SOURCE CODE, BUT WITH future.apply TO FIX PARALLEL COMPUTATION
    # https://github.com/xzhoulab/iDEA/blob/a30fc8ec655813678049c5b4b4b0aa6a21c944f0/R/iDEASummary.R#L103
    # here setting default parameters from tutorial at: https://xzhoulab.github.io/iDEA/
    object = idea
    init_tau = c(-2,0.5)
    em_iter = 15
    mcmc_iter = 1000
    min_degene = 0
    #### start snippet from iDEA source
    num_gene <- object@num_gene
    num_annot <- length(object@annotation)
    init_beta <- object@summary[,1]
    # this line is the only deviation from the iDEA sourcecode;
    # we use the future.apply package because the original code, using 'pbmclapply'
    # 1) occasionally hangs in WSL and 2) doesn't perform well on Windows
    # (simple workaround as we only use this specific iDEA configuration)
    # (note that we can't print to console from within future "workers", so iDEA-specific-errors are lost / not reported)
    res_idea <- future.apply::future_sapply(1:num_annot, FUN = function(x) {
      set.seed(x)
      Annot <- rep(0, object@num_gene)
      Annot[object@annotation[[x]]] <- 1
      Annot <- Annot - mean(Annot)
      Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
      t1 <- system.time(model1 <- try( res <- iDEA::EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
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
    idea = iDEA::iDEA.louis(idea)

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
    # If the geneset algorithm returned a non-finite pvalue, we consider it "tested but found not signif"
    # (note that leaving NA values as-is would affect subsequent multiple testing correction, which in turn would be a
    # source of bias if many  "no effect whatsoever" genesets are returned as NA pvalues by the algorithm)
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
