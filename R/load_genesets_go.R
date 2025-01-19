
#' Load GO annotations via Bioconductor packages (e.g. org.Hs.eg.db for Human)
#'
#' @description
#' Download and import genesets from the GO database using the Bioconductor infrastructure.
#' Use the `goat::load_genesets_go_fromfile` function for more fine-grained control over the GO database version that you use; it allows you to import NCBI gene2go files
#'
#' @details
#' Note that org.Hs.eg.db pulls data semi-annually from NCBI gene2go (as do similar databases for other species),
#' but the GO database version returned by this function is tied to the version of the org.Hs.eg.db on your computer (this is controlled by the Bioconductor infrastructure).
#' E.g. in an outdated R installation you may get outdated GO annotations as well.
#'
#' The actual GO database version that is retrieved is returned by this function in the `source_version` column.
#' @param include_child_annotations boolean; include annotations against child terms? In most situations, TRUE (default) is the desired setting
#' @param taxid integer value that indicates the taxonomy id. Default: 9606 (Human, Homo sapiens). Currently supported options:
#'
#' * 9606 = Human (Homo sapiens)
#' * 9598 = Chimpanzee (Pan troglodytes)
#' * 9544 = Rhesus monkey (Macaca mulatta)
#' * 10090 = Mouse (Mus musculus)
#' * 10116 = Rat (Rattus norvegicus)
#' * 7955 = Zebrafish (Danio rerio)
#' * 7227 = Fruit fly (Drosophila melanogaster)
#' * 6239 = Worm (Caenorhabditis elegans)
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_go_bioconductor = function(include_child_annotations = TRUE, taxid = 9606) {
  genes = go_id = GOID = TERM = ONTOLOGY = go_domain = source_version = go_name = ngenes = parent_id = child_id = relation = name = NULL # fix invisible bindings R package NOTE

  ### input validation
  stopifnot("include_child_annotations parameter should be a single boolean value" = length(include_child_annotations) == 1 && include_child_annotations %in% c(TRUE, FALSE))
  stopifnot("taxid parameter should be a single positive integer value" = length(taxid) == 1 && is.numeric(taxid) && is.finite(taxid) && taxid > 0)

  taxid_lookup = taxonomy_identifiers()
  taxid = as.integer(taxid) # enforce integer type, also deals with rounding etc.
  stopifnot("provided taxid parameter is not supported by this function, c.f. function parameter documentation" = taxid %in% taxid_lookup$taxonomy_id)

  taxid_info = taxid_lookup[match(taxid, taxid_lookup$taxonomy_id), ] # always yields exactly 1 row, because upstream we validated taxid is in this table
  stopifnot("failed to match input taxid against the table from taxonomy_identifiers()" = nrow(taxid_info) == 1)
  check_dependency("AnnotationDbi", "loading GO data via Bioconductor")
  check_dependency("GO.db", "loading GO data via Bioconductor")
  check_dependency(taxid_info$bioconductor_package, paste("loading", taxid_info$common_name, "GO species database via Bioconductor"))

  # get the database object from the respective bioconductor package
  ORG_ES_DB = eval(parse(text = paste0(taxid_info$bioconductor_package, "::", taxid_info$bioconductor_package))) # this is basically; ORG_ES_DB=org.Hs.eg.db::org.Hs.eg.db
  # dev note: alternatively, use require() then get() but note that this may fail in some cases, e.g. after updating packages in the same session

  ### Bioconductor GO terms and respective annotations (Entrez gene format)

  # using Bioconductor package 'GO.db', get a table of GO ID, name, ontology (CC/BP/MF)
  go_terms = suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), columns = c("TERM","ONTOLOGY"), keytype = "GOID", multiVals = "first"))
  # using Bioconductor package 'org.Hs.eg.db', get a list of GO ID with values entrez ID
  keys = setdiff(unique(go_terms$GOID), c("GO:0003674", "GO:0008150", "GO:0005575")) # exclude top-level ontologies like "molecular function" and CC/BP counterparts (basically entire realm)
  # bugfix; previously we used `AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "GO")` to extract all unique GO term IDs but this seems to be bugged; it leaves out many legit terms (e.g. ribosomal subunit GO:0044391) that do have annotations in org.Hs.eg.db::org.Hs.eg.db !
  # extract annotations
  go_annotations_entrez = suppressMessages(AnnotationDbi::mapIds(ORG_ES_DB, keys = keys, column = "ENTREZID", keytype = ifelse(include_child_annotations, "GOALL", "GO"), multiVals = "list"))
  # extract all metadata from the org.XX.eg.db package
  go_annotations_metadata_table = AnnotationDbi::metadata(ORG_ES_DB)
  go_annotations_metadata_table = go_annotations_metadata_table[
    !is.na(go_annotations_metadata_table$name) & nchar(go_annotations_metadata_table$name) >= 2 &
      !is.na(go_annotations_metadata_table$value) & nchar(go_annotations_metadata_table$value) >= 2, ] # remove rows without data/information
  # find the organism in the source data. Not strictly needed because we know the taxonomy, but good validation if there are problems with source data
  # if available use ORGANISM, fallback is SPECIES, otherwise enter the requested species
  source_organism = stats::na.omit(go_annotations_metadata_table$value[match(c("ORGANISM", "SPECIES"), toupper(go_annotations_metadata_table$name))])
  if(length(source_organism) == 0) {
    source_organism = taxid_info$scientific_name
  }
  # in older versions of all species databases, GOSOURCEDATE is always present but unfortunately more recent versions do not include the actual GO source date !
  # For example, org.Mm.eg.db version 2024-Sep20 has an empty value for the GOSOURCEDATE field !
  source_version = stats::na.omit(go_annotations_metadata_table$value[match(c("GOSOURCEDATE", "GOEGSOURCEDATE", "EGSOURCEDATE", "UPSOURCEDATE"), toupper(go_annotations_metadata_table$name))])
  if(length(source_version) == 0) {
    source_version = "version/date not described in Bioconductor package"
    print(go_annotations_metadata_table)
  }

  go_annotations_metadata = paste0(taxid_info$bioconductor_package, " - ", source_organism[1], " - ", source_version[1]) # take first element from all available organism/version data


  ### convert Bioconductor data into a table compatible with this R package

  result = tibble::tibble(go_id = rep(names(go_annotations_entrez), lengths(go_annotations_entrez)),
                          genes = unlist(go_annotations_entrez, recursive = FALSE, use.names = FALSE)) |>
    # enforce entrez gene IDs to be integers by stripping non-numeric parts
    # (not strictly needed atm, just a safeguard against future upstream changes, e.g. prefixing entrez IDs with 'entrez:')
    mutate(genes = gsub("\\D+","", genes)) |>
    filter(genes != "") |>
    mutate(genes = as.integer(genes)) |>
    # remove duplicate goterm*genes entries, if any
    distinct(go_id, genes, .keep_all = TRUE) |>
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
    rename(id = go_id, name = go_name)


  ### ontology DAG

  extract_links = function(GODBLINKS, relation_accept) {
    # extract direct parents
    GOdata = as.list(GODBLINKS)
    # unlist the named vector per GO term into a long-format table
    links = dplyr::bind_rows(sapply(names(GOdata), function(n)
      data.frame(child_id=n, parent_id=unname(GOdata[[n]]), relation=names(GOdata[[n]]), row.names = NULL), simplify = FALSE, USE.NAMES = FALSE))
    links |>
      # remove unsupported relation types
      filter(relation %in% relation_accept) |>
      # retain only unique parent/child links, discarding relation types
      distinct(parent_id, child_id) |>
      tibble::as_tibble()
  }

  # accepted relation types (as specified in GO.obo) + variations that the Bioconductor GO.db might use (e.g. "isa" or "part of")
  relation_accept = c("is_a","part_of", "regulates", "positively_regulates", "negatively_regulates")
  relation_accept = unique(c(relation_accept, sub("_", "", relation_accept), sub("_", " ", relation_accept)))
  # links from GO.db::GOCCPARENTS are supposed to be within-GO-domain (e.g. no links from BP to CC)
  links_cc = extract_links(GO.db::GOCCPARENTS, relation_accept)
  links_bp = extract_links(GO.db::GOBPPARENTS, relation_accept)
  links_mf = extract_links(GO.db::GOMFPARENTS, relation_accept)


  ### compose final result

  result = result |>
    left_join(
      bind_rows(links_cc, links_bp, links_mf) |>  tidyr::chop(cols = parent_id),
      by = c("id"="child_id")
    ) |>
    # column ordering and renaming
    select(source, source_version, id, name, parent_id, genes, ngenes)

  attr(result, "settings") <- sprintf("load_genesets_go_bioconductor(include_child_annotations=%s) = '%s'",
                                      include_child_annotations, go_annotations_metadata)

  message(paste("load_genesets_go_bioconductor(): data version =", go_annotations_metadata))
  return(result)
}



#' construct a geneset table from gene2go and OBO files
#'
#' @description
#' This function is used to load Gene Ontology (GO) genesets from files that you
#' manually downloaded from the links below. This enables the use of the latest data
#' from GO (in contrast,  Bioconductor GO data may lag behind current data considerably).
#' To construct genesets from available raw data, download the "gene2go" file
#' (the gene annotations) from below NCBI link and download the GO OBO
#' (ontology terms and relations to respective parent/child terms)  from below
#' geneontology.org link. Provide the full path to the downloaded file to this function.
#' Both "gzipped" and "uncompressed" files are supported.
#'
#' We encourage you to rename the files after your downloaded them such that
#' the date of download in incorporated; this ensures you can always keep track of
#' the GO database version that was used! For example, rename the downloaded
#' "gene2go.gz" file to "gene2go_2024-01-31.gz".
#'
#' Download link for gene2go file; https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#'
#' Download link for gene ontology OBO file; http://current.geneontology.org/ontology/go.obo
#'
#' @examples
#'   # TODO: update the filenames to your downloaded files
#'   file_gene2go = "C:/DATA/gene2go_2024-01-01.gz"
#'   file_goobo = "C:/DATA/go_2024-01-01.obo"
#'   if(file.exists(file_gene2go) && file.exists(file_goobo)) {
#'     genesets_asis = load_genesets_go_fromfile(file_gene2go, file_goobo)
#'   }
#' @param file_gene2go full path to the gene2go file from NCBI. Also works with a gzipped file; gene2go.gz
#' @param file_goobo full path to the OBO file from geneontology.org. Also works with a gzipped file; obo.gz
#' @param include_child_annotations boolean; include annotations against child terms? In most situations, TRUE (default) is the desired setting
#' @param taxid_filter return gene annotations only for the specified taxonomy identifier (integer value). Typical options are 9606 (Human, the default) or 10090 (mouse). Importantly, select the taxonomy/organism that is also in your input gene list (which typically contains human Entrez gene identifiers)
#' @return table with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
#' @export
load_genesets_go_fromfile = function(file_gene2go, file_goobo, include_child_annotations = TRUE, taxid_filter = 9606) {
  genes = parent_id_recursive = namespace = source_version = name = ngenes = parent_id = NULL # fix invisible bindings R package NOTE

  # input validation
  stopifnot("include_child_annotations parameter should be a single boolean value" = length(include_child_annotations) == 1 && include_child_annotations %in% c(TRUE, FALSE))
  stopifnot("taxid_filter parameter should be a single positive integer value" = length(taxid_filter) == 1 && is.numeric(taxid_filter) && is.finite(taxid_filter) && taxid_filter > 0)
  stopifnot("file_gene2go parameter should be a single string (full path to respective file)" = length(file_gene2go) == 1 && !is.na(file_gene2go) && is.character(file_gene2go) && nchar(file_gene2go) > 0)
  stopifnot("file_goobo parameter should be a single string (full path to respective file)" = length(file_goobo) == 1 && !is.na(file_goobo) && is.character(file_goobo) && nchar(file_goobo) > 0)
  if(!file.exists(file_gene2go)) stop(paste("cannot find input file:", file_gene2go))
  if(!file.exists(file_goobo)) stop(paste("cannot find input file:", file_goobo))

  # parse input data
  gene2go = go_gene2go(file_gene2go, taxid_filter = taxid_filter)
  go = go_obo(file_goobo, rename_namespace = TRUE, remove_obsolete = TRUE)

  # gene*term pairs (direct annotations)
  annotations = gene2go |> select(id, genes) |> tidyr::unchop(genes)

  if(include_child_annotations) {
    # long-format table of all parent terms
    id_to_parents = go |>
      select(id, parent_id_recursive) |>
      tidyr::unchop(parent_id_recursive) # drop empty is fine; ignore terms without parents

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

    # debug;
    # id_to_parents |> filter(parent_id_recursive == "GO:0005615") |> left_join(go |> select(id, name, namespace, obsolete)) |> arrange(name) |> print(n=Inf)
    # id_to_parents |> filter(parent_id_recursive == "GO:0005576") |> left_join(go |> select(id, name, namespace, obsolete)) |> arrange(name) |> print(n=Inf)
    # annotations |> filter(id %in% c("GO:0005576", "GO:0005615") & genes == 21) # this gene needs to be rolled up from 0005615 to 0005576
    # annotations_rollup |> filter(genes == 21 & id %in% c("GO:0005576", "GO:0005615"))
  }

  # prepare a result table in the same format as other geneset import functions
  result = go |>
    mutate(
      source = paste0("GO_", namespace),
      source_version = paste(file_gene2go, file_goobo)
    ) |>
    select(source, source_version, id, name, parent_id) |>
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
#' @description note that this file lacks parent/child relations, so we only learn 'direct annotations'
#' @param f full path to gene2go file stored on the computer, e.g. previously downloaded from https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#' @param taxid_filter return gene annotations only for the specified taxonomy identifier (integer value). Typical options are 9606 (Human, the default) or 10090 (mouse). Importantly, select the taxonomy/organism that is also in your input gene list (which typically contains human Entrez gene identifiers)
#' @return a tibble with columns; source, source_version, id, name, genes, ngenes
#' @export
go_gene2go = function(f, taxid_filter = 9606) {
  taxid = qualifier = goid = goterm = geneid = category = genes = source_version = name = ngenes = NULL # fix invisible bindings R package NOTE
  stopifnot("taxid_filter parameter should be a single positive integer value" = length(taxid_filter) == 1 && is.numeric(taxid_filter) && is.finite(taxid_filter) && taxid_filter > 0)
  stopifnot("f parameter should be a single string (full path to respective file)" = length(f) == 1 && !is.na(f) && is.character(f) && nchar(f) > 0)
  if(!file.exists(f)) stop(paste("cannot find input file:", f))

  gene2go = data.table::fread(f, data.table = FALSE, showProgress = FALSE)
  # for robustness against future changes, e.g. colnames changing 'case'
  colnames(gene2go) = gsub("[^a-z]+", "", tolower(colnames(gene2go)))

  # only taxid-of-interest
  gene2go = gene2go |> filter(taxid == taxid_filter)

  gene2go |>
    # remove negative annotation. e.g. "NOT part_of" (do this after taxid filter for efficiency)
    filter( ! grepl("^not ", qualifier, ignore.case = TRUE)) |>
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
#' @param child array of GO term IDs that represent children
#' @param parent array of GO term IDs that represent respective parents
#' @noRd
go_find_parents = function(child, parent) {
  newchild = NULL # fix invisible bindings R package NOTE
  stopifnot(length(child) > 0 & length(child) == length(parent) & !anyNA(child) & !anyNA(parent))

  # initialize result matrix; provided child/parent pairs (direct links)
  input = data.frame(child, parent)
  x = data.frame(child, parent, done = ! parent %in% child)

  max_depth = 1000
  n = 0
  while(any(!x$done)) {
    n = n + 1
    if(n > max_depth) {
      stop("unexpected depth while traversing GO structure, is the provided structure not a DAG or tree ?!")
    }

    # construct map for TODO rows. Use left-join so we include 1:n links (i.e. don't use match())
    y = left_join(data.frame(child = x$parent[!x$done], newchild = x$child[!x$done]), input, by = "child") |>
      select(child = newchild, parent)
    # only done if new parent/end was already followed through
    y$done = ! y$parent %in% x$child

    # flag previous batch as done
    x$done = TRUE
    # update resultset
    x = bind_rows(x, y)
  }

  # remove duplicates and return as a tibble with unique child/parent pairs
  x$done = NULL
  tibble::as_tibble(x) |> distinct_all()
}



#' simple vectorized parsing of GO OBO file without any dependencies (beyond dplyr/tibble/tidyr)
#'
#' @description note that we remove links between GO terms that are across GO domains (e.g. no CC to MF relations)
#' The only supported relations are those that match this regex;
#' `"^(is_a:|relationship: part_of|relationship: regulates|relationship: positively_regulates|relationship: negatively_regulates)"`
#'
#' @param f full path to go.obo file stored on the computer, e.g. previously downloaded from http://current.geneontology.org/ontology/go.obo . Also works with a gzipped file; obo.gz
#' @param rename_namespace boolean; rename official namespace values like 'cellular_component' to CC? (analogous for BP and MF)
#' @param remove_obsolete boolean; remove obsoleted terms?
#' @return tibble with ontology terms and their relations
#' @export
go_obo = function(f, rename_namespace = TRUE, remove_obsolete = TRUE) {
  group = lines = isterm = isid = isname = isdef = isnamespace = isobsolete = isrelationship = isparentlink = NULL # fix invisible bindings R package NOTE
  id = name = definition = namespace = obsolete = qc = parent_id = parent_namespace = parent = parent_id_recursive = NULL # fix invisible bindings R package NOTE

  stopifnot("rename_namespace parameter should be a single boolean value" = length(rename_namespace) == 1 && rename_namespace %in% c(TRUE, FALSE))
  stopifnot("remove_obsolete parameter should be a single boolean value" = length(remove_obsolete) == 1 && remove_obsolete %in% c(TRUE, FALSE))
  stopifnot("f parameter should be a single string (full path to respective file)" = length(f) == 1 && !is.na(f) && is.character(f) && nchar(f) > 0)
  if(!file.exists(f)) stop(paste("cannot find input file:", f))

  go = tibble::tibble(
    # instead of `readLines(f, warn = FALSE)` we here use data.table to read all lines (`sep=NULL` to read everything in a single column)
    # this is fast, has the advantage that we can read .gz compressed input files, does not require additional dependencies (we already use data.table in other functions)
    lines = data.table::fread(f, sep = NULL, colClasses = "character", header = FALSE, encoding = "UTF-8", stringsAsFactors = FALSE, data.table = FALSE)[,1],
    # lines = readLines(f, warn = FALSE),
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
  golinks_recursive = go_find_parents(child = golinks$id, parent = golinks$parent_id)

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
