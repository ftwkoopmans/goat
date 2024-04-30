
#' Parse HGNC gene identifier lookup table that was downloaded from genenames.org into a table with HGNC ID, symbol, synonym (NA if unavailable), entrez ID
#'
#' download link: https://www.genenames.org/download/statistics-and-files/
#' table: "Complete dataset download links" -->> "Complete HGNC approved dataset" -->> download the "TXT" table
#' filename is typically something like hgnc_complete_set.txt
#' URL; https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
#'
#' # alternatively;
#' table: "Total Approved Symbols" -->> "TXT" / "text file in TSV format"
#' filename is typically something like non_alt_loci_set.txt
#'
#' @param filename full path to the downloaded table (expected to be tsv format, typically has .txt or .tsv file extension)
#' @return a long-format table with columns; hgnc_id, hgnc_symbol, type, value
#' @export
hgnc_idmap_table = function(filename) {
  symbol = hgnc_id = hgnc_symbol = entrez_id = prev_symbol = alias_symbol = NULL # fix invisible bindings R package NOTE
  stopifnot("parameter filename must be an existing file" = length(filename) == 1 && is.character(filename) && file.exists(filename))

  # parse HGNC table
  hgnc = data.table::fread(filename, data.table = FALSE, stringsAsFactors = FALSE)
  colnames(hgnc) = tolower(colnames(hgnc))
  cols_expect = c("hgnc_id", "symbol", "alias_symbol", "prev_symbol", "entrez_id")
  cols_alt = c("HGNC ID", "Approved symbol", "Previous symbols", "Alias symbols", "NCBI Gene ID")
  ok = all(cols_expect %in% colnames(hgnc))
  if(ok) {
    hgnc = hgnc[,cols_expect]
  } else {
    ok = all(cols_alt %in% colnames(hgnc))
    if(!ok) stop(paste("expected a HGNC table with columns;", paste(cols_expect, collapse = ", ")))
    hgnc = hgnc[,match(cols_alt, colnames(hgnc))]
    colnames(hgnc) = cols_expect
  }

  hgnc = hgnc |>
    rename(hgnc_symbol = symbol) |>
    # remove invalid rows; no id, no valid symbol, no entrez ID
    filter(!is.na(hgnc_id) & hgnc_id != "" & !is.na(hgnc_symbol) & nchar(hgnc_symbol) >= 2 & !is.na(entrez_id) & grepl("^\\d+", as.character(entrez_id))) |>
    mutate(
      hgnc_symbol = toupper(hgnc_symbol),
      prev_symbol = ifelse(is.na(prev_symbol), "", prev_symbol),
      alias_symbol = ifelse(is.na(alias_symbol), "", alias_symbol),
      synonyms = paste(prev_symbol, alias_symbol),
      entrez_id = as.integer(gsub("\\D.*", "", as.character(entrez_id))) # remove multiple, if any
    )

  # parse synonyms into a long-format lookup table & remove empty / ambiguous
  # (synonyms that overlap with main symbol OR are found 2+ times)
  l = strsplit(toupper(hgnc$synonyms), "[ ,;|]+") # support various delimiters
  map_synonym = data.frame(
    hgnc_id = rep(hgnc$hgnc_id, lengths(l)),
    symbol = unlist(l, recursive = FALSE, use.names = FALSE)
  ) |>
    filter(nchar(symbol) >= 2 &  ! symbol %in% hgnc$hgnc_symbol) |>
    # deal with duplicate synonyms within 1 row/hgnc_id
    # (e.g. due to our case conversion or overlap between previous and alias columns)
    distinct_all()

  # count how often each synonym is found --> remove those with multiple entries
  symbol_ambiguous = map_synonym |> count(symbol) |> filter(n > 1) |> pull(symbol)
  map_synonym = map_synonym |> filter( ! symbol %in% symbol_ambiguous)

  hgnc |> select(hgnc_id, hgnc_symbol, entrez_id) |> left_join(map_synonym |> rename(synonym = symbol), by = "hgnc_id")
}



#' Map the the symbol column in a table to HGNC human gene IDs by matching official gene symbols and synonyms
#'
#' @examples
#'   # TODO: update the filename to your downloaded file
#'   # download instructions in the documentation of `hgnc_idmap_table()`
#'   f = "C:/DATA/HGNC/hgnc_complete_set.txt"
#'
#'   if(file.exists(f)) {
#'     df = data.frame(symbol = c("vamp2", "STXBP1", "UNC18", NA, "PSD95", "NOT-A-GENE"))
#'     hgnc = hgnc_idmap_table(f)
#'     df = symbol_to_entrez(df, hgnc)
#'     print(df)
#'   }
#' @param x a data.table with a column symbol
#' @param hgnc HGNC lookup table from `hgnc_idmap_table()`
#' @return entrez gene IDs are returned in the "gene" column of table `x`. Additionally, columns "entrez_id", "hgnc_id" and "hgnc_symbol"
#' @export
symbol_to_entrez = function(x, hgnc) {
  symbol = NULL # fix invisible bindings R package NOTE
  if(!is.data.frame(x) || nrow(x) == 0 || !"symbol" %in% colnames(x) || !is.character(x$symbol) ) {
    stop("x must be a data.frame with column 'symbol' (character type)")
  }
  if(!is.data.frame(hgnc) || nrow(hgnc) == 0 || !all(c("hgnc_id", "hgnc_symbol", "synonym", "entrez_id") %in% colnames(hgnc)) ) {
    stop("hgnc must be a data.frame with columns 'hgnc_id', 'hgnc_symbol', 'synonym', 'entrez_id' as typically prepared using the hgnc_idmap_table() function")
  }

  x$symbol_input = x$symbol
  rows_symbol_fail = is.na(x$symbol) | nchar(x$symbol) < 2
  x$symbol = toupper(x$symbol)

  # match directly by official gene symbol
  i = match(x$symbol, hgnc$hgnc_symbol)
  x$hgnc_id = hgnc$hgnc_id[i]

  # for rows that failed to match, try matching synonyms
  rows = is.na(x$hgnc_id)
  if(any(rows)) {
    i = match(x$symbol[rows], hgnc$synonym)
    x$hgnc_id[rows] = hgnc$hgnc_id[i]
  }

  # erase whatever happened to invalid input
  x$hgnc_id[rows_symbol_fail] = NA

  # add hgnc_symbol and entrez ID
  i = match(x$hgnc_id, hgnc$hgnc_id)
  x$hgnc_symbol = hgnc$hgnc_symbol[i]
  x$entrez_id = hgnc$entrez_id[i]
  x$gene = x$entrez_id

  # restore input symbols
  x$symbol = x$symbol_input
  x$symbol_input = NULL

  # log status to console. use distinct() because some input tables might contain multiple entries
  tmp = x |> filter(rows_symbol_fail == FALSE) |> distinct(symbol, .keep_all = TRUE)
  n_input = nrow(tmp)
  n_fail = sum(is.na(tmp$hgnc_id))
  message(sprintf("%d / %d (%.1f%%) unique symbols could not be mapped to a HGNC human gene ID",
              n_fail, n_input, n_fail / n_input * 100 ))

  return(x)
}
