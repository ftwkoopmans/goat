
#' validate that a geneset table is compatible with this R package
#'
#' @param genesets geneset tibble to validate, e.g. results from `load_genesets_go_bioconductor()`
#' @param require_signif should we check for "ngenes_signif" column ?
#' @param check_gene_identiers optionally skip the validation of gene identifiers
validate_genesets = function(genesets, require_signif = TRUE, check_gene_identiers = TRUE) {
  # 1) data.frame with all required columns
  ok = is.data.frame(genesets) &&
    nrow(genesets) > 0 &&
    all(c("source", "source_version", "id", "name", "genes", "ngenes") %in% colnames(genesets)) &&
    (!require_signif || "ngenes_signif" %in% colnames(genesets) )


  # 2) check column types
  if(ok) {
    types = sapply(genesets, typeof)
    ok = all(c("source", "source_version", "id", "name", "genes", "ngenes") %in% names(types)) &&
      all(types[c("source", "source_version", "id", "name")] == "character") &&
      types["genes"] == "list" &&
      types["ngenes"] == "integer"

    if(require_signif) {
      ok = ok && "ngenes_signif" %in% names(types) && types["ngenes_signif"] == "integer"
    }
  }

  if(!ok) {
    err = "geneset table should be a non-empty data.frame/tibble with these columns (and types); source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)"
    if(require_signif) {
      err = paste0(err, "ngenes_signif (int)\nDid you forget to use the filter_genesets() prior ?")
    }
    stop(err)
  }


  # 3) check non-list columns for NA values
  ok = !anyNA(genesets) &&
    all(genesets$source != "") &&
    all(genesets$source_version != "") &&
    all(genesets$id != "") &&
    all(genesets$name != "") &&
    all(is.finite(genesets$ngenes) & is.integer(genesets$ngenes)) &&
    (!require_signif || all(is.finite(genesets$ngenes_signif) & is.integer(genesets$ngenes_signif)) )


  if(!ok) {
    stop("geneset table should not contain empty strings or NA/infinite values")
  }


  # 4) validate gene values; have to be non-empty characters
  if(check_gene_identiers) {
    x_aslist = utils::as.relistable(genesets$genes)
    # ignore warning about recursive unlisting; we validate below to enforce 'no recursive lists allowed'
    x = unlist(genesets$genes, recursive = F, use.names = F)
    if(length(x) > 0) {
      ok = FALSE
      if(is.character(x)) {
        ok = !anyNA(x) && all(x != "")
      }
      if(is.integer(x) || is.numeric(x)) { # also allow 'numeric' to relax compatability a bit
        ok = all(is.finite(x))
      }
      if(!ok) {
        stop("geneset table should contain a 'genes' column that is of list-type and contains only characters (not NA, not '') or only integer gene IDs (not NA, only finite values)")
      }

      # convert numeric non-integer gene identifiers, code documentation/rationale @ validate_genelist()
      if( ! is.integer(x) && is.numeric(x)) {
        gene_as_int = as.integer(round(x, digits = 0))
        if(any(abs(x - gene_as_int) > 10^-6)) {
          stop("geneset table should contain whole numbers (integers) in the 'genes' column, if numeric identifiers are provided")
        }

        # note that we have to map the unlisted genes back into the same list structure
        tmp = utils::relist(gene_as_int, skeleton = x_aslist)
        class(tmp) = "list" # remove the 'relistable' class / typedef
        genesets$genes = tmp
      }
    }
  }

  # TODO: gene identifier validation for genes_signif, if present, analogous to above

  return(genesets)
}
