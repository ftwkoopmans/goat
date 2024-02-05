
#' Plot a treemap
#'
#' @description simple wrapper around the treemap R package. To customize this plot, copy/paste its code and tweak parameters as desired
#' @examples \dontrun{
#' # use GOAT as per usual
#' genelist = goat_example_datasets[["Hondius 2021:mass-spec:PMID33492460"]]
#' genesets = load_genesets_go_bioconductor()
#' genesets_filtered = filter_genesets(genesets, genelist)
#' goat_result = test_genesets(
#'   genesets_filtered, genelist, method = "goat", score_type = "effectsize",
#'   padj_method = "bonferroni", padj_cutoff = 0.05
#' )
#' # subset GO CC results
#' x = goat_result |> filter(signif & source == "GO_CC")
#' tm = treemap_data(
#'   geneset_ids = x$id,
#'   genesets = genesets,
#'   genesets_test_result = x,
#'   simplify = "leaf_only" # options: none/leaf_only/prune_singletons/pvalue
#' )
#' treemap_plot(tm$treemap_plotdata)
#' }
#' @param x `treemap_plotdata` data table that was computed by the `treemap_data()` function
#' @param label_group set TRUE to show only group-level labels
#' @export
treemap_plot = function(x, label_group = FALSE) {
  check_dependency("treemap", "plotting treemaps")
  if(!is.data.frame(x) || nrow(x) == 0) {
    return()
  }
  # use ggplot color palette, skipping the last few colors (too close to initial red hue)
  nclr = n_distinct(x$group_name)
  clr = rev(utils::head(gg_color_hue(nclr + ceiling(nclr / 10)), n = nclr))
  # add alpha to mute the colors
  if(nclr > 8) {
    clr = paste0(clr, utils::head(rep(c("BB", "99", "77"), each = ceiling(nclr/3)), n = nclr))
  } else {
    clr = paste0(clr, "88")
  }
  # debug, plot color palette; barplot(rep(1, nclr), col = clr)

  treemap::treemap(
    x,
    index = c("group_name", "subgroup_name"),
    vSize = "ngenes_sqrt",
    title = "",
    title.legend = "",
    position.legend = "none",
    fontsize.labels = c(ifelse(label_group, 10, 0), ifelse(label_group, 0, 10)),
    fontcolor.labels = "#000000",
    bg.labels = 240,
    lowerbound.cex.labels = 0.1,
    vColor = "group_name",
    type = "categorical",
    palette = clr
    # palette = "HCL",
    # palette.HCL.options = list(hue_start=30, hue_end=360, luminance = 90)
  )
}



#' Construct tree and treemap data structures from geneset parent/child relations
#'
#' @examples \dontrun{
#'  # refer to the goat::treemap_plot() function for a complete example
#' }
#' @param geneset_ids vector of geneset identifiers
#' @param genesets entire geneset table; typically the complete GO database
#' @param genesets_test_result geneset testing results; the output from `test_genesets()`
#' @param simplify strategy for reducing the genesets returned in the treemap. Options;
#' "leaf_only" (most stringent, returns only leafs in the tree structure)
#' "prune_singletons" (remove parent terms that have exactly 1 child)
#' "pvalue" (remove parent terms where the child term p-value is at least 4 times better)
#' "none" (default; return all significant genesets that are not a "grouping term" in the treemap)
#' @param toplevel_max_ngenes groups in the treemap should not have more than this many genes ('ngenes' in geneset test results). If not set, this defaults to 50% of the total number of unique genes in the geneset test results
#' @export
treemap_data = function(geneset_ids, genesets, genesets_test_result, simplify = "none", toplevel_max_ngenes = NA) {
  ods = ontology_data_structures(geneset_ids = geneset_ids, genesets = genesets, genesets_test_result = genesets_test_result)
  if(length(toplevel_max_ngenes) == 1 && is.na(toplevel_max_ngenes)) {
    ngenes_total = n_distinct(unlist(genesets_test_result$genes))
    toplevel_max_ngenes = ngenes_total * 0.5
  }
  build_treemap(ods = ods, simplify = simplify, toplevel_max_ngenes = toplevel_max_ngenes)
}



#' depth-first recursion
#'
#' @param p p
#' @param edgelist edgelist
#' @param result result
edgelist_find_children_recursive = function(p, edgelist, result) {
  children = setdiff(edgelist$child[edgelist$parent == p], result)
  if(length(children) == 0) return(result)
  for(cld in children) {
    result = c(result, cld, edgelist_find_children_recursive(cld, edgelist, result))
  }
  return(result)
}



#' Within each recursive path; halt at first colored node
#'
#' @param x geneset ID
#' @param geneset_table data.frame with columns 'id' and 'parent_id'
#' @param query set of all colored nodes (vector)
#' @param result all colored parents, accumulated recursively (vector)
edgelist_find_colored_parents_recursive = function(x, geneset_table, query, result) {
  i = match(x, geneset_table$id)
  if(!is.na(i)) {
    parents = geneset_table$parent_id[[i]]
    if(length(parents) > 0) { # at least 1 parent
      parents_colored = intersect(parents, query) # can be NULL
      parents_rest = setdiff(parents, parents_colored) # can be NULL
      # add colored parents to results
      result = c(result, parents_colored)
      # recurse for the others
      for(p in parents_rest) {
        result = edgelist_find_colored_parents_recursive(p, geneset_table, query, result)
      }
    }
  }
  return(result)
}



#' Convert edgelist table into a nested list structure
#'
#' @param x element ID
#' @param name element name
#' @param ngenes element ngenes property
#' @param edgelist data.frame with properties parent, child, child_name, child_ngenes
edgelist_as_nested = function(x, name, ngenes, edgelist) {
  l = list(id = x, name = name, ngenes = ngenes, children = list())
  rows = edgelist$parent == x
  if(any(rows)) {
    el = edgelist[!rows,]
    for(i in which(rows)) {
      l$children[[length(l$children) + 1]] = edgelist_as_nested(edgelist$child[i], edgelist$child_name[i], edgelist$child_ngenes[i], el)
    }
  }
  return(l)
}



#' For each element in a DAG, count the number of recursive child terms and those that are on shortlist
#'
#' @param x list representing a nested DAG
#' @param shortlist vector of element IDs
nested_aggregate_child_stats = function(x, shortlist) {
  x$recursiveChildren = c()
  x$recursiveChildrenShortlist = c()
  x$recursiveChildren_length = x$recursiveChildrenShortlist_length = 0
  if(length(x$children) > 0) {
    for(i in seq_along(x$children)) {
      x$children[[i]] = nested_aggregate_child_stats(x$children[[i]], shortlist)
      x$recursiveChildren = c(x$recursiveChildren, x$children[[i]]$id, x$children[[i]]$recursiveChildren)
    }
    x$recursiveChildren        = unique(x$recursiveChildren) # dedupe; DAGs may yield duplicates
    x$recursiveChildren_length = length(x$recursiveChildren)
    x$recursiveChildrenShortlist        = intersect(x$recursiveChildren, shortlist)
    x$recursiveChildrenShortlist_length = length(x$recursiveChildrenShortlist)
  }
  return(x)
}



#' Recursively replace a DAG element by its children while the child has the same gene count as the parent
#'
#' @param x list representing a nested DAG
nested_find_equivalent_child_recursive = function(x) {
  if(length(x$children) > 0) {
    child_ngenes = unlist(lapply(x$children, "[[", "ngenes"))
    # iterate in order of least number of ngenes (most specific term first = sort in ascending order)
    for(i in order(child_ngenes, decreasing = FALSE)) {
      # found a replacement
      if(x$children[[i]]$recursiveChildrenShortlist_length == x$recursiveChildrenShortlist_length) {
        return(nested_find_equivalent_child_recursive(x$children[[i]]))
      }
    }
  }
  return(x)
}



#' Traverse DAG and recursively split parent term into child terms until none are larger than N genes
#'
#' @param obj list representing a nested DAG
#' @param threshold stop if `obj$ngenes <= threshold`
#' @param result list of resulting elements (when calling this function, use default empty list)
nested_find_level1_children = function(obj, threshold, result = list()) {
  # halt recursion when there are no children, or node has less than <threshold> genes
  if(length(obj$children) == 0 || obj$ngenes <= threshold) {
    result[[length(result) + 1]] = obj
    return(result)
  }
  for(child in obj$children) {
    result = nested_find_level1_children(child, threshold, result)
  }
  return(result)
}



#' Construct DAG between genesets-of-interest as edgelist
#'
#' @description this function reduces the complete DAG to a subset with only parameter genesets
#' @param ids shortlist of geneset IDs
#' @param genesets importantly, genesets should be the input genesets and not the "filtered" genesets because
#' only the former contains the complete ontological structure (parent/child links between genesets)
edgelist_from_ontology = function(ids, genesets) {
  if(length(ids) == 0) {
    return()
  }
  edgelist = NULL
  for(q in ids) { # q = "GO:0097060"
    x = edgelist_find_colored_parents_recursive(q, genesets, ids, NULL)
    if(length(x) == 0) {
      x = "root" # default root
    }
    edgelist = bind_rows(edgelist, tibble::tibble(child = q, parent = x))
  }
  edgelist |> distinct_all()
}



#' Find root ID and name in edgelist
#'
#' throw error if there is no single root
#'
#' @param edgelist result from e.g. `edgelist_from_ontology()`
#' @param genesets filtered genesets
edgelist_find_root = function(edgelist, genesets) {
  root_id = setdiff(unique(edgelist$parent), unique(edgelist$child))
  stopifnot(length(root_id) == 1)
  i = match(root_id, genesets$id)
  root_name = ifelse(!is.na(genesets$name[i]), genesets$name[i], "root")
  root_ngenes = ifelse(!is.na(genesets$ngenes[i]), genesets$ngenes[i], 0)
  return(list(id = root_id, name = root_name, ngenes = root_ngenes))
}



#' Compute tree and treemap data structures from geneset tables that include parent/child links
#'
#' @param geneset_ids vector of geneset identifiers
#' @param genesets entire geneset table; typically the complete GO database
#' @param genesets_test_result geneset testing results; the output from `test_genesets()`
ontology_data_structures = function(geneset_ids, genesets, genesets_test_result) {
  name = ngenes = pvalue = parent_ngenes = child = NULL # fix invisible bindings R package NOTE
  stopifnot("parameter geneset_ids must not be empty and contain valid genesets_test_result$id values" = length(geneset_ids) > 0 && is.character(geneset_ids) && all(geneset_ids %in% genesets_test_result$id))
  stopifnot("parameter genesets must be a geneset table. This function requires ontological structure in the genesets, represented by a 'parent_id' column (e.g. does NOT work with genesets imported in GMT format)" =
              is.data.frame(genesets) && nrow(genesets) > 0 && "parent_id" %in% colnames(genesets))
  stopifnot("parameter genesets must be a geneset table. This function requires ontological structure in the genesets, represented by a 'parent_id' column (e.g. does NOT work with genesets imported in GMT format)" =
              is.data.frame(genesets_test_result) && nrow(genesets_test_result) > 0 && "parent_id" %in% colnames(genesets_test_result))

  # DAG: edgelist = from parent-child links in genesets-of-interest to a connected edgelist ('between colored nodes')
  dag_edgelist = edgelist_from_ontology(geneset_ids, genesets)
  dag_edgelist_root = edgelist_find_root(dag_edgelist, genesets_test_result)
  # add metadata from geneset testing results (ngenes filtered + pvalue)
  dag_edgelist = dag_edgelist |>
    left_join(genesets_test_result |> select(child = id, child_name = name, child_ngenes = ngenes, child_pvalue = pvalue), by = "child") |>
    left_join(genesets_test_result |> select(parent = id, parent_name = name, parent_ngenes = ngenes, parent_pvalue = pvalue), by = "parent")

  # DAG: nested = from edgelist to nested data structure (incl. recursive child term counts etc.)
  dag_nested = edgelist_as_nested(dag_edgelist_root$id, dag_edgelist_root$name, dag_edgelist_root$ngenes, dag_edgelist)

  # tree: edgelist = reduce the DAG edgelist to only 'best' link per child (sort ascending by parent_ngenes)
  tree_edgelist = dag_edgelist |>
    arrange(parent_ngenes) |>
    distinct(child, .keep_all = TRUE)

  # tree: nested = analogous to DAG (same root even), but with reduced edgelist as input
  tree_nested = edgelist_as_nested(dag_edgelist_root$id, dag_edgelist_root$name, dag_edgelist_root$ngenes, tree_edgelist)
  # TODO: can we sort by geneset*geneset similarity at each (or top 3~5) level of the tree?

  return(list(geneset_ids = geneset_ids, dag_edgelist=dag_edgelist, dag_nested=dag_nested, tree_edgelist=tree_edgelist, tree_nested=tree_nested))
}



#' Compute treemap data structures
#'
#' @description
#' - shortlist = genesets/DAG-nodes to plot
#' - update nested DAG data structure with recursive counts of shortlist elements
#' - find level-1 elements to start with (possibly further down than direct children of root)
#' - collapse children of level-1 elements (either from tree structure or from greedy aggregation starting at largest element)
#' - restructure treemap data into format suitable for plotting
#' @param ods result from `ontology_data_structures()`
#' @param simplify strategy for reducing the genesets returned in the treemap. Options;
#' "leaf_only" (most stringent, returns only leafs in the tree structure)
#' "prune_singletons" (remove parent terms that have exactly 1 child)
#' "pvalue" (remove parent terms where the child term p-value is at least 4 times better)
#' "none" (default; return all significant genesets that are not a "grouping term" in the treemap)
#' @param toplevel_max_ngenes groups in the treemap should not have more than this many genes ('ngenes' in geneset test results)
build_treemap = function(ods, simplify = "none", toplevel_max_ngenes = Inf) {
  parent = parent_pvalue = child = child_pvalue = child_name = child_ngenes = children = name = ngenes = group_ngenes = subgroup_name = NULL # fix invisible bindings R package NOTE
  stopifnot(length(simplify) == 1 && simplify %in% c("leaf_only", "prune_singletons", "pvalue", "none"))
  stopifnot(length(toplevel_max_ngenes) == 1 && is.numeric(toplevel_max_ngenes) && !is.na(toplevel_max_ngenes) && toplevel_max_ngenes >= 0)


  ## shortlist = genesets/DAG-nodes to plot
  shortlist = ods$geneset_ids
  if(simplify == "leaf_only") {
    shortlist = setdiff(unique(ods$dag_edgelist$child), unique(ods$dag_edgelist$parent))
  } else if(simplify == "prune_singletons") {
    # only retain nodes that have more or less (leafs) than 1 child
    shortlist = setdiff(ods$geneset_ids, ods$dag_edgelist |> count(parent) |> filter(n==1) |> pull(parent))
  } else if(simplify == "pvalue") {
    # remove parents where pvalues is not at least X times better than child pvalue
    factor_parent_must_be_better = 4
    shortlist = setdiff(
      ods$geneset_ids,
      ods$dag_edgelist |> filter(is.finite(child_pvalue) & is.finite(parent_pvalue) & child_pvalue * factor_parent_must_be_better <= parent_pvalue) |> pull(parent)
    )
  }

  ## update nested DAG data structure with recursive counts of shortlist elements
  dag = nested_aggregate_child_stats(ods$dag_nested, shortlist)

  ## find level-1 elements to start with (possibly further down than direct children of root)
  # start with the children of the root; collect nearest children that have fewer than <threshold> genes
  obj_toplevel = list()
  for(child_level1 in dag$children) {
    for(tmp in nested_find_level1_children(obj = child_level1, threshold = toplevel_max_ngenes)) {
      obj_toplevel[[length(obj_toplevel) + 1]] = tmp
    }
  }
  # take unique set and sort by largest-first
  obj_toplevel = obj_toplevel[!duplicated(unlist(lapply(obj_toplevel, "[[", "id")))]
  obj_toplevel = obj_toplevel[order(unlist(lapply(obj_toplevel, "[[", "ngenes")), decreasing = TRUE)]

  ## collapse children of level-1 elements (either from tree structure or from greedy aggregation starting at largest element)
  treemap_data = list()
  shortlist_done = NULL
  # greedy aggregation, start with level-1 term that has most genes
  for(obj in obj_toplevel) {
    i_shortlist = i_root = NULL
    # current term already covered
    if(obj$id %in% shortlist_done) {
      next
    }
    # current term is on the shortlist of IDs we're looking for
    if(obj$id %in% shortlist) {
      i_root = obj
      i_shortlist = i_root$id
      # include child terms that are on the shortlist, if any
      if(i_root$recursiveChildrenShortlist_length > 0) {
        i_shortlist = c(i_shortlist, i_root$recursiveChildrenShortlist)
      }
      i_shortlist = setdiff(i_shortlist, shortlist_done) # skip already done
    } else {
      # current term is not on the shortlist
      # find most-specific child that covers the same number of leafs
      i_shortlist = setdiff(obj$recursiveChildrenShortlist, shortlist_done) # skip already done
      if(length(i_shortlist) > 0) {
        # deal with top-level terms that have only 1 child (which may also have only 1 child, etc.)
        i_root = nested_find_equivalent_child_recursive(obj)
      }
    }

    # update results
    if(length(i_shortlist) > 0) {
      treemap_data[[length(treemap_data) + 1]] = list(id = i_root$id, name = i_root$name, ngenes = i_root$ngenes, children = i_shortlist)
      shortlist_done = c(shortlist_done, i_shortlist)
    }
  }


  ### convert the nested result data into a table suitable for treemap plotting
  treemap_plotdata = NULL
  if(length(treemap_data) > 0) {
    treemap_plotdata = bind_rows(treemap_data) |>
      rename(group = id, subgroup = children, group_name = name, group_ngenes = ngenes) |>
      left_join(ods$dag_edgelist |> select(subgroup = child, subgroup_name = child_name, ngenes = child_ngenes), by = "subgroup") |>
      mutate(ngenes_sqrt = sqrt(ngenes)) |>
      # (current iteration of) above code already reduces the DAG to a tree.
      # To ensure future compatibility we here enforce child uniqueness again
      arrange(desc(group_ngenes)) |>
      distinct(subgroup_name, .keep_all = TRUE)
  }

  cat(nrow(treemap_plotdata), "/", length(ods$geneset_ids), " genesets remain after simplifying ontology structure by '", simplify , "'\n", sep="")
  return(list(treemap_data = treemap_data, treemap_plotdata = treemap_plotdata))
}
