
#' parse genesets from the SynGO database
#'
#' Workflow;
#' - obtain the input file from; https://www.syngoportal.org
#' - click "bulk download SynGO release ..." for SynGO release of interest
#' - unzip
#' - call this function with the full file path to the 'syngo_ontologies.xlsx' file
#' @examples
#'   # TODO: update the filename to your downloaded file
#'   f = "C:/DATA/SynGO_bulk_download_release_20210225/syngo_ontologies.xlsx"
#'   if(file.exists(f)) genesets_asis = load_genesets_syngo(f)
#' @param filename full path to the "syngo_ontologies.xlsx" file that was extracted from a SynGO bulk download ZIP archive
#' @param gene_database gene IDs to return. must be any of; "entrez" (default), "hgnc", "ensembl"
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_syngo = function(filename, gene_database = "entrez") {
  name = domain = genes = go_name = go_id = go_domain = source_version = ngenes = parent_id = NULL # fix invisible bindings R package NOTE
  stopifnot("filename parameter must be an existing file" = length(filename) == 1 && is.character(filename) && file.exists(filename))
  stopifnot("gene_database parameter must be any of; 'entrez' or 'hgnc' or 'ensembl'" = length(gene_database) == 1 && gene_database %in% c("hgnc", "entrez", "ensembl"))

  col_gene = "hgnc_id"
  if(gene_database == "entrez") col_gene = "entrez_id"
  if(gene_database == "ensembl") col_gene = "ensembl_id"

  input = readxl::read_excel(filename)
  # robust colnames for compat with v1.0
  colnames(input) = gsub("[^a-z0_]+", "", tolower(colnames(input)))
  colnames(input)[colnames(input) == "gotermid"] = "id"
  colnames(input)[colnames(input) == "gotermname"] = "name"
  colnames(input)[colnames(input) == "godomain"] = "domain"
  colnames(input)[colnames(input) == "goparenttermid"] = "parent_id"
  colnames(input)[colnames(input) == "geneshgnc_id"] = "hgnc_id"
  # colnames(input)[colnames(input) == "geneshgnc_symbol"] = "hgnc_symbol"

  cols_required = c("id", "name", "domain", "parent_id", col_gene)
  if(!all(cols_required %in% colnames(input))) {
    stop(paste("missing column names:", paste(setdiff(cols_required, colnames(input)), collapse = " , ")))
  }


  result = input |>
    select(go_id = id, go_name = name, go_domain = domain, parent_id, genes = {{col_gene}}) |>
    filter(nchar(genes) > 0) |>
    mutate( # remove ontology ID from name column
      go_name = sub(" *\\((SYNGO|GO):.*", "", go_name),
      # split delimited genes into a list and unnest the column (so we can do vectorized mutate and filter)
      genes = strsplit(genes, " *[,;] *")
    ) |>
    # to long format
    tidyr::unnest(genes)

  if(gene_database == "hgnc") {
    result = result |> # remove non-integer characters from ID
      mutate(genes = paste0("HGNC:", gsub("\\D+", "", genes))) |>
      filter(genes != "HGNC:") # no empty IDs !
  }
  # HGNC/entrez; remove invalid = retain only integer IDs
  if(gene_database == "entrez") {
    result = result |> filter(grepl("^\\d+$", genes)) |> mutate(genes = as.integer(genes))
  }
  # ensembl; simply check for string length, trust otherwise
  if(gene_database == "ensembl") {
    result = result |> filter(nchar(genes) > 5)
  }

  result = result |> # no duplicate goterm*gene entries (there should be none in SynGO, but enforce in code anyway)
    distinct(go_id, genes, .keep_all = TRUE) |>
    arrange(genes) |>
    tidyr::chop(genes) |>
    mutate(source = paste0("SYNGO_", go_domain),
           source_version = filename,
           ngenes = lengths(genes)) |>
    select(source, source_version, id=go_id, name=go_name, parent_id, genes, ngenes) # column ordering and renaming

  attr(result, "settings") <- sprintf("load_genesets_syngo(filename='%s', gene_database='%s')",
                                      filename, gene_database)
  return(result)
}
