
#' parse genesets in GMT format where gene identifiers are numeric Entrez gene IDs
#'
#' # Example data;
#' URL: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C5
#' download this data: KEGG subset of curated pathways -->> NCBI (Entrez) Gene IDs
#' filename should be something like "c2.cp.kegg.v2023.1.Hs.entrez.gmt"
#'
#' @examples
#'   # TODO: update the filename to your downloaded file
#'   f = "C:/DATA/c2.cp.kegg.v2023.1.Hs.entrez.gmt"
#'   if(file.exists(f)) genesets_asis = load_genesets_gmtfile(f, label = "KEGG")
#' @param filename input file for this function should be the full path to genesets defined in GMT format
#' @param label a shortname for the genesets in this file, for example "GO_CC", "KEGG", "MY_DB_V1". This will be stored in the 'source' column of the resulting table. Importantly, multiple testing correction in GOAT is grouped by this 'source' column so you probably want to use a different label for each collection-of-genesets that you load. Must not be empty, only allowed characters are; upper/lower-case letter, numbers 0-9 and underscore
#' @return tibble with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_gmtfile = function(filename, label) {
  stopifnot("filename parameter must be an existing file" = length(filename) == 1 && is.character(filename) && file.exists(filename))
  stopifnot("label parameter must be at least 2 characters (and 100 at most) and may only contain alphanumeric characters and underscore (e.g. no space, '-', semicolons, etc.)" =
              length(label) == 1 && nchar(label) >= 2 && nchar(label) <= 100 && grepl("^[a-zA-Z0-9_]+$", label))

  ### expected format:
  # reach line is a geneset, separated by tabs
  # the first element to be the geneset ID, second element is optionally some description/string, rest are gene identifiers
  x = strsplit(gsub("(^ +)|( +$)", "", readLines(filename)), " *\t *") # split by tab while also trimming of excess whitespace (if any)
  # remove lines with less than 3 elements
  x = x[lengths(x) >= 3]
  # extract ID; first element
  x_id = unlist(lapply(x, utils::head, n = 1))
  # x_name = unlist(lapply(x, "[", 2)); x_name[is.na(x_name) | x_name == ""] = x_id[is.na(x_name) | x_name == ""]
  # retain only genes = purge the first 2 elements
  x = lapply(x, utils::tail, n = -2)
  # specific for GOAT pipeline, assume all integer IDs
  x = lapply(x, as.integer)
  # enforce unique (probably not needed, but no harm in double-checking)
  x = lapply(x, unique)
  # store in a table that follows the format expected in this package
  result = tibble::tibble(source = label, source_version = filename, id = x_id, name = x_id,
                          # important technical detail; cast to a typed list for consistency with other load_genesets_... functions
                          genes = vctrs::as_list_of(x), ngenes = lengths(x)) |>
    # QC double-check; if there are multiple entries for the same geneset ID, just keep the first
    distinct(id, .keep_all = TRUE)

  attr(result, "settings") <- sprintf("load_genesets_gmtfile(filename='%s', label='%s')", filename, label)
  return(result)
}
