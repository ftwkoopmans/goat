
#' human gene (NCBI entrez ID) annotations from the GO database using the 'org.Hs.eg.db' Bioconductor package
#'
#' note that org.Hs.eg.db pulls data semi-annually from NCBI gene2go,
#' the GO database version returned by this function is tied to the version of the org.Hs.eg.db on your computer.
#' Note that the actual GO database version is returned in the result column.
#'
#' @examples \dontrun{
#'   genesets_asis = load_genesets_go_bioconductor()
#' }
#' @param include_child_annotations boolean; include annotations against child terms? In most situations, TRUE (default) is the desired setting
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_go_bioconductor = function(include_child_annotations = TRUE) {
  genes = go_id = GOID = TERM = ONTOLOGY = go_domain = source_version = go_name = ngenes = NULL # fix invisible bindings R package NOTE
  check_dependency("AnnotationDbi", "loading GO data via Bioconductor")
  check_dependency("GO.db", "loading GO data via Bioconductor")
  check_dependency("org.Hs.eg.db", "loading GO data via Bioconductor")

  # using Bioconductor package 'GO.db', get a table of GO ID, name, ontology (CC/BP/MF)
  go_terms = suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), columns = c("TERM","ONTOLOGY"), keytype = "GOID", multiVals = "first"))
  # using Bioconductor package 'org.Hs.eg.db', get a list of GO ID with values entrez ID
  keys = setdiff(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "GO"), c("GO:0003674", "GO:0008150", "GO:0005575")) # exclude top-level ontologies like "molecular function" and CC/BP counterparts (basically entire realm)
  go_annotations_entrez = suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = keys, column = "ENTREZID", keytype = ifelse(include_child_annotations, "GOALL", "GO"), multiVals = "list"))
  go_annotations_entrez = go_annotations_entrez[lengths(go_annotations_entrez) > 0]
  # GO DB version
  go_annotations_metadata = AnnotationDbi::metadata(org.Hs.eg.db::org.Hs.eg.db)
  go_annotations_metadata = paste(go_annotations_metadata$value[match(c("GOSOURCENAME", "ORGANISM", "GOSOURCEDATE"), go_annotations_metadata$name)], collapse = " - ")

  # from geneset list to long-format table
  result = tibble::tibble(go_id = rep(names(go_annotations_entrez), lengths(go_annotations_entrez)),
         genes = unlist(go_annotations_entrez, recursive = F, use.names = F)) |>
    # enforce entrez gene IDs to be integers by stripping non-numeric parts
    # (not strictly needed atm, just a safeguard against future upstream changes, e.g. prefixing entrez IDs with 'entrez:')
    mutate(genes = gsub("\\D+","", genes)) |>
    filter(genes != "") |>
    mutate(genes = as.integer(genes)) |>
    # remove duplicate goterm*genes entries, if any
    distinct(go_id, genes, .keep_all = T) |>
    # back to list format
    tidyr::chop(cols = genes) |>
    # add goterm metadata; name and domain
    left_join(tibble::as_tibble(go_terms) |> select(go_id = GOID, go_name = TERM, go_domain = ONTOLOGY), by = "go_id") |>
    # package-specific metadata we'll use downstream
    mutate(
      source = paste0("GO_", go_domain),
      source_version = go_annotations_metadata,
      ngenes = lengths(genes)
    ) |>
    select(source, source_version, id = go_id, name = go_name, genes, ngenes) # column ordering and renaming


  attr(result, "settings") <- sprintf("load_genesets_go_bioconductor(include_child_annotations=%s) = '%s'",
                                      include_child_annotations, go_annotations_metadata)
  return(result)
}



#' construct a geneset table from gene2go and OBO files
#'
#' Download link for gene2go file; https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#' Download link for gene ontology OBO file; http://current.geneontology.org/ontology/go.obo
#'
#' @examples \dontrun{
#'   genesets_asis = load_genesets_go_fromfile(
#'     file_gene2go = "C:/DATA/download_2023-01-01/gene2go",
#'     file_goobo = "C:/DATA/download_2023-01-01/go.obo"
#'   )
#' }
#' @param file_gene2go gene2go file from NCBI. Also works with the gzipped file gene2go.gz
#' @param file_goobo OBO file from geneontology.org
#' @param include_child_annotations boolean; include annotations against child terms? In most situations, TRUE (default) is the desired setting
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_go_fromfile = function(file_gene2go, file_goobo, include_child_annotations = TRUE) {
  genes = parent_id_recursive = namespace = source_version = name = ngenes = NULL # fix invisible bindings R package NOTE
  # parse input data
  gene2go = go_gene2go(file_gene2go, taxid_filter = 9606)
  go = go_obo(file_goobo, rename_namespace = TRUE, remove_obsolete = TRUE)

  # gene*term pairs (direct annotations)
  annotations = gene2go |> select(id, genes) |> tidyr::unchop(genes)

  if(include_child_annotations) {
    # long-format table of all parent terms
    id_to_parents = go |>
      select(id, parent_id_recursive) |>
      tidyr::unchop(parent_id_recursive)

    # gene*term pairs for all 'hierarchical rollup', associating genes with all recursive parents
    annotations_rollup = id_to_parents |>
      # remove terms without annotations to speedup the 'join'
      filter(id %in% annotations$id) |>
      # left_join by annotated term --> now a gene annotated against term X is also associated with all parent terms of X
      left_join(annotations, by = "id")

    # finally, combine the direct and 'rollup' annotations and retain only the unique set
    # note that this overwrites the current 'direct annotations' table
    annotations = bind_rows(annotations, annotations_rollup |> select(id = parent_id_recursive, genes)) |>
      distinct(id, genes)
  }

  # prepare a result table in the same format as other geneset import functions
  result = go |>
    mutate(
      source = paste0("GO_", namespace),
      source_version = paste(file_gene2go, file_goobo),
    ) |>
    select(source, source_version, id, name) |>
    # note; drop_na() not strictly needed, nor is the sorting of genes
    left_join(annotations |> tidyr::drop_na() |> arrange(genes) |> tidyr::chop(cols = genes), by = "id") |>
    mutate(ngenes = lengths(genes)) |>
    # don't return empty goterms (no annotations)
    filter(ngenes > 0)

  rm(gene2go, go, annotations, id_to_parents, annotations_rollup)
  gc(full = TRUE, verbose = FALSE)


  attr(result, "settings") <- sprintf("load_genesets_go_fromfile(file_gene2go='%s', file_goobo='%s', include_child_annotations='%s')",
                                      file_gene2go, file_goobo, include_child_annotations)
  return(result)
}



#' parse gene2go file
#'
#' note that it lacks parent/child relations, so from this file we only learn 'direct annotations'
#'
#' @param f full path to gene2go file stored on the computer, e.g. previously downloaded from https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#' @param taxid_filter taxonomy id, integer
go_gene2go = function(f, taxid_filter = 9606) {
  taxid = qualifier = goid = goterm = geneid = category = genes = source_version = name = ngenes = NULL # fix invisible bindings R package NOTE

  gene2go = data.table::fread(f, data.table = FALSE, showProgress = FALSE)
  # for robustness against future changes, e.g. colnames changing 'case'
  colnames(gene2go) = gsub("[^a-z]+", "", tolower(colnames(gene2go)))

  gene2go |>
    # only taxid-of-interest
    filter(taxid == taxid_filter) |>
    # remove negative annotation. e.g. "NOT part_of" (do this after taxid filter for efficiency)
    filter( ! grepl("^not ", qualifier, ignore.case = T)) |>
    as_tibble() |>
    # retain only data we need & rename columns
    select(id = goid, name = goterm, genes = geneid, source = category) |>
    # enforce unique
    distinct(id, genes, .keep_all = TRUE) |>
    tidyr::chop(cols = genes) |>
    mutate(
      # for robustness against future changes, match against lower-case and try multiple strings
      source = tolower(source),
      source = paste0(
        "GO_",
        ifelse(source %in% c("component", "cc"), "CC",
               ifelse(source %in% c("process", "bp"), "BP",
                      ifelse(source %in% c("function", "mf"), "MF", NA_character_)))
      ),
      source_version = f,
      ngenes = lengths(genes)
    ) |>
    tidyr::drop_na() |>
    select(source, source_version, id, name, genes, ngenes) # column ordering
}



#' from a directed edgelist (child-to-parent) to a recursive lookup of the final parent term / root
#'
#' instead of recursive lookups per term, we iteratively perform vectorized matching over the entire table
#' since the GO tree/DAG is quite shallow, this is very fast
#'
#' @param child array of GO term IDs that represent children
#' @param parent array of GO term IDs that represent respective parents
go_find_parents = function(child, parent) {
  stopifnot(length(child) > 0 & !anyNA(child) & !anyNA(parent) & length(child) == length(parent))

  # initialize result matrix; provided child/parent pairs
  mat = cbind(child, parent)
  colnames(mat) = c("child", "parent")
  last = parent

  max_depth = 1000
  n = 1
  while(TRUE) {
    n = n + 1
    if(n > max_depth) {
      stop("unexpected depth while traversing GO structure, is the provided structure not a DAG or tree ?!")
    }

    # match GO ID from 'latest parent list' against 'child', so we know what the parent-of-parent is
    j = match(last, child)
    if(any(!is.na(j))) {
      # new result = link between 'latest parent list' and their respective parents
      tmp = cbind(child, parent[j])
      # remove NA matches (e.g. previous parent was NA to begin with, or there is no parent-of-parent)
      tmp = tmp[!is.na(tmp[,1]) & !is.na(tmp[,2]), , drop=F]
      # keep going on success, break otherwise
      if(nrow(tmp) > 0) {
        mat = rbind(mat, tmp)
        last = parent[j]
      } else {
        break
      }
    } else {
      break
    }
  }

  # remove duplicates
  # our vectorized matching is quite fast but generates duplicates because we pass the same 'path to tree root' many times
  mat = mat[!duplicated(mat), , drop=F]

  # return as a tibble
  tibble::as_tibble(mat)
}



#' simple vectorized parsing of GO OBO file without any dependencies (beyond dplyr/tibble/tidyr)
#'
#' note that we remove links between GO terms that are across GO domains (e.g. no CC to MF relations)
#'
#' @param f full path to go.obo file stored on the computer, e.g. previously downloaded from http://current.geneontology.org/ontology/go.obo
#' @param rename_namespace boolean; rename official namespace values like 'cellular_component' to CC? (analogous for BP and MF)
#' @param remove_obsolete boolean; remove obsoleted terms?
go_obo = function(f, rename_namespace = TRUE, remove_obsolete = TRUE) {
  group = lines = isterm = isid = isname = isdef = isnamespace = isobsolete = isrelationship = isparentlink = NULL # fix invisible bindings R package NOTE
  id = name = definition = namespace = obsolete = qc = parent_id = parent_namespace = parent = parent_id_recursive = NULL # fix invisible bindings R package NOTE

  go = tibble::tibble(
    lines = readLines(f, warn = F),
    group = cumsum(as.integer(grepl("^\\[", lines))),
    isterm = lines == "[Term]",
    isid = grepl("^id:", lines),
    isname = grepl("^name:", lines),
    isdef = grepl("^def:", lines),
    isnamespace = grepl("^namespace:", lines),
    isobsolete = lines == "is_obsolete: true",
    # isisa = grepl("^is_a:", lines),
    isrelationship = grepl("^relationship:", lines),
    # all link fields; c("is_a", "regulates", "part_of", "has_part", "happens_during", "negatively_regulates", "positively_regulates", "occurs_in", "ends_during")
    # not included; occurs_in / happens_during / ends_during / has_part (careful, this is a parent->child link, all others are child->parent)
    isparentlink = grepl("^(is_a:|relationship: part_of|relationship: regulates|relationship: positively_regulates|relationship: negatively_regulates)", lines)
  ) |>
    # for each 'group of lines', enforce that the group describes exactly 1 [Term] header, 1 id, 1 name and 1 namespace
    group_by(group) |>
    mutate(qc = sum(isterm) == 1 & sum(isid) == 1 & sum(isname) == 1 & sum(isnamespace)) |>
    ungroup() |>
    filter(qc == TRUE)


  # compose basic GO term info table
  goterms = go |>
    group_by(group) |>
    summarise(
      id = lines[isid],
      name = lines[isname],
      definition = lines[isdef],
      namespace = lines[isnamespace],
      obsolete = any(isobsolete),
      .groups = "drop"
    ) |>
    mutate(
      id = sub("^id: *", "", id),
      name = sub("^name: *", "", name),
      definition = sub("^def: *", "", definition),
      definition = gsub('(^"|"$|" *\\[.*)', "", definition),
      namespace = sub("^namespace: *", "", namespace)
    ) |>
    # remove root term and obsolete terms (optional)
    filter( ! name %in% c("cellular_component", "biological_process", "molecular_function") & (!remove_obsolete | !obsolete) )


  if(rename_namespace) {
    goterms$namespace[goterms$namespace == "cellular_component"] = "CC"
    goterms$namespace[goterms$namespace == "biological_process"] = "BP"
    goterms$namespace[goterms$namespace == "molecular_function"] = "MF"
  }

  # for each GO term, collect relations
  golinks = go |>
    filter(isparentlink) |>
    select(lines, group) |>
    # fast regex to grab the first GO ID;
    # remove everything up to the first number, then reconstruct 'GO:' prefix
    mutate(parent_id = sub("^(\\d+).*", "GO:\\1", gsub("^\\D+", "", lines))) |>
    # should not occur, but remove all lines where we failed to find a proper GO ID
    # also remove go terms that were removed upstream (e.g. root terms or obsoleted terms)
    filter(parent_id != lines & group %in% goterms$group & parent_id %in% goterms$id) |>
    distinct(group, parent_id) |>
    # collect GO Term info for the current line = child term
    left_join(goterms |> select(group, id, namespace), by = "group") |>
    # collect GO Term info for the identified link to parent
    left_join(goterms |> select(id, parent_namespace = namespace), by = c("parent_id"="id")) |>
    # remove links across domains, e.g. MF term where 'parent link' leads to a CC term
    filter(!is.na(namespace) & !is.na(parent_namespace) & parent_namespace == namespace) |>
    select(id, parent_id)

  # find recursive path up to root / last parent
  golinks_recursive = go_find_parents(golinks$id, golinks$parent_id)

  ### QC
  # # show all OBO lines related to term 'synaptic vesicle'
  # go |> filter(group %in% (goterms |> filter(name == "synaptic vesicle") |> pull(group) ))
  # # show all parent terms for term 'synaptic vesicle'
  # golinks_recursive |> filter(child %in% (goterms |> filter(name == "synaptic vesicle") |> pull(id) )) |> left_join(goterms |> select(id, name), by = c(parent="id"))


  goterms |>
    select(-group) |>
    left_join(golinks |> tidyr::chop(cols = parent_id), by = "id") |>
    left_join(golinks_recursive |> rename(parent_id_recursive = parent) |> tidyr::chop(cols = parent_id_recursive), by = c("id"="child"))
}
