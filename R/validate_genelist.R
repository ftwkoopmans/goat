
#' validate a table with genes (that should be tested in overrepresentation-analysis) for compatibility with this R package
#'
#' @param genelist gene tibble to validate
validate_genelist = function(genelist) {
  # 1) data.frame with all required columns
  ok = is.data.frame(genelist) &&
    nrow(genelist) > 0 &&
    all(c("gene", "signif") %in% colnames(genelist))

  # 2) check column types
  if(ok) {
    types = sapply(genelist, typeof)
    ok = all(c("gene", "signif") %in% names(types)) &&
      types["gene"] %in% c("character", "integer", "numeric", "double") &&
      types["signif"] == "logical"
  }
  if(!ok) {
    stop("genelist table should be a data.frame/tibble with these columns (and types); gene (character or integer), signif (logical/boolean)")
  }

  # 3) check for NA or empty-string values  (note that the 'signif' column can be NA)
  ok = FALSE
  if(is.character(genelist$gene)) {
    ok = !anyNA(genelist$gene) && all(genelist$gene != "")
  }
  if(is.integer(genelist$gene) || is.numeric(genelist$gene)) { # also allow 'numeric' to relax compatability a bit
    ok = all(is.finite(genelist$gene)) # disallow NA and Inf
  }
  if(!ok) {
    stop("genelist table should not contain empty/missing/NA/Inf values in the 'gene' column")
  }

  # 4) for numeric type gene identifiers we allow double-type but require integer values --> validate & type conversion
  # this is just for user convenience, ideally we'd require strict value types but that might force users into
  # type conversion mistakes so we'll handle this here.
  # e.g. user might have integer IDs from NCBI Entrez, but unwittingly store them in a numeric/double column type
  # here we test that conversion to integer doesn't lose information (i.e. numeric-type IDs have no decimal values)
  if( ! is.integer(genelist$gene) && is.numeric(genelist$gene)) {
    # before conversion to int, we should round() because integer conversion drops all decimals (rounds down)
    # see also the help/documentation @ ?as.integer
    # testcase;  tmp = 1/(1-0.99); sprintf("%.14f", tmp); as.integer(tmp); all.equal(as.integer(tmp), tmp)
    gene_as_int = as.integer(round(genelist$gene, digits = 0))
    if(any(abs(genelist$gene - gene_as_int) > 10^-6)) { # test with some minor tolerance for imprecision
      stop("genelist table should contain whole numbers (integers) in the 'gene' column, if numeric identifiers are provided")
    }
    # alternatively, check without tolerance;
    # stopifnot(genelist$gene %% 1 == 0); genelist$gene = as.integer(genelist$gene)

    # type conversion
    genelist$gene = gene_as_int
  }

  # 5) genes cannot be duplicated
  if(anyDuplicated(genelist$gene)) {
    stop("genelist table should not contain duplicate values in the 'gene' column")
  }

  return(genelist)
}
