
#' parse genesets from the SynGO database
#'
#' Workflow;
#' - obtain the input file from; https://www.syngoportal.org
#' - click "bulk download SynGO release ..." for SynGO release of interest
#' - unzip
#' - call this function with the full file path to the 'syngo_ontologies.xlsx' file
#' @examples \dontrun{
#'   genesets_asis = load_genesets_syngo(
#'     "C:/DATA/SynGO_bulk_download_release_20210225/syngo_ontologies.xlsx"
#'   )
#' }
#' @param filename input file for this function should be the full path name for "syngo_ontologies.xlsx" from a SynGO bulk download
#' @param gene_database gene IDs to return. must be any of; "entrez" (default), "hgnc", "ensembl"
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_syngo = function(filename, gene_database = "entrez") {
  name = domain = genes = go_name = go_id = go_domain = source_version = ngenes = NULL # fix invisible bindings R package NOTE
  stopifnot("filename parameter must be an existing file" = length(filename) == 1 && is.character(filename) && file.exists(filename))
  stopifnot("gene_database parameter must be any of; 'entrez' or 'hgnc' or 'ensembl'" = length(gene_database) == 1 && gene_database %in% c("hgnc", "entrez", "ensembl"))

  col_gene = "hgnc_id"
  if(gene_database == "entrez") col_gene = "entrez_id"
  if(gene_database == "ensembl") col_gene = "ensembl_id"

  input = readxl::read_excel(filename)
  stopifnot(c("id", "name", "domain", col_gene) %in% colnames(input))

  result = input |>
    select(go_id = id, go_name = name, go_domain = domain, genes = {{col_gene}}) |>
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
    distinct(go_id, genes, .keep_all = T) |>
    arrange(genes) |>
    tidyr::chop(genes) |>
    mutate(source = paste0("SYNGO_", go_domain),
           source_version = filename,
           ngenes = lengths(genes)) |>
    select(source, source_version, id=go_id, name=go_name, genes, ngenes) # column ordering and renaming

  attr(result, "settings") <- sprintf("load_genesets_syngo(filename='%s', gene_database='%s')",
                                      filename, gene_database)
  return(result)
}
