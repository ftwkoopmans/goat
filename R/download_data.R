
#' Download the datasets that were used in the GOAT manuscript
#'
#' @description
#' Downloads OMICs-based datasets that were used in the GOAT manuscript from the GOAT GitHub page.
#' This file is cached in the output directory and only needs to be downloaded once. Multiple datasets
#' are included and their names include the respective PubMed identifiers (PMID).
#'
#' If you encounter technical difficulties, try to;
#'
#' 1) download the file by copy/pasting this URL into your browser: https://github.com/ftwkoopmans/goat/raw/main/analyses/goat_manuscript_datasets.rda
#' 2) load the data in R using the following 2 lines of code, here assuming you stored the downloaded file at C:/data/goat_manuscript_datasets.rda
#'
#' `load("C:/data/goat_manuscript_datasets.rda")`
#'
#' `genelist = goat_manuscript_datasets.rda[["Wingo 2020:mass-spec:PMID32424284"]]`
#'
#' @param output_dir full path to the directory where the downloaded files should be stored. Directory is created if it does not exist.
#' e.g. `output_dir="~/data"` on unix systems, `output_dir="C:/data"` on Windows, or set to `output_dir=getwd()` to write output to the current working directory
#' @param ignore_cache boolean, set to TRUE to force re-download and ignore cached data, if any. Default: FALSE
#' @return a list of genelist data tables. The names of the list represent the datasets,
#' values in the list are data tables that can be used as a "genelist" in the GOAT R package
#' @export
download_goat_manuscript_data = function(output_dir, ignore_cache = FALSE) {
  stopifnot("parameter output_dir must be a single string and represent a directory on your computer" = length(output_dir) == 1 && is.character(output_dir) && !is.na(output_dir))
  stopifnot("parameter ignore_cache must be a single boolean value" = length(ignore_cache) == 1 && ignore_cache %in% c(TRUE, FALSE))
  sprintf_template_downloadfail = "failed to download %s and store it at %s\nTry an alternative output_dir parameter or follow the download_goat_manuscript_data() function documentation to manually download the file and load it in R"

  # create dir if it does not exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if(!dir.exists(output_dir)) {
      stop(paste0("Could not create the requested output directory: ", output_dir, "\nTry to provide an existing directory as parameter for 'output_dir'"))
    }
  }

  # load from cache, or download if not available
  filename = paste0(output_dir, "/goat_manuscript_datasets.rda")
  url = "https://github.com/ftwkoopmans/goat/raw/main/analyses/goat_manuscript_datasets.rda"
  if(ignore_cache || !file.exists(filename)) {
    message(paste("downloading", url, "..."))
    utils::download.file(url, filename, mode = "wb")
    if(!file.exists(filename)) {
      stop(sprintf(sprintf_template_downloadfail, url, filename))
    }
    message(paste0("downloaded data was stored at: ", filename))
  } else {
    message(paste0("cached data was retrieved from: ", filename))
  }

  # load RData file into environment
  e = new.env()
  load(filename, envir = e)

  # validate that the expected variable is present
  if(!is.list(e$goat_manuscript_datasets)) {
    stop("failed to load RData file; it did not contain expected variable 'goat_manuscript_datasets'")
  }

  return(e$goat_manuscript_datasets)
}



#' Download and parse geneset collections from the GOAT GitHub repository
#'
#' @description while the Bioconductor respository is extensive, contains data for many species and is a part of
#' a larger infrastructure, it might contain outdated GO data when the user is not using the latest R version.
#' If users are on an R version that is a few years old, so will the GO data from Bioconductor be.
#'
#' As an alternative, we store gene2go data from NCBI (for Human genes only!) at the GOAT GitHub repository.
#' This function allows for a convenient way to download this data and then parse the genesets.
#'
#' Alternatively you can browse the files in the data branch of the GOAT GitHub repository and download these files manually,
#' then load them via the GOAT R function `load_genesets_go_fromfile()`.
#' To view these data, open this URL in a browser; https://github.com/ftwkoopmans/goat/tree/data
#' You can also use this R package to see all available data via `available_genesets_goatrepo()`
#'
#' By default (empty `version` parameter), this function will first check the online GOAT GitHub repository to find the
#' most recent version/date, then download the respective data.
#'
#' @examples \donttest{
#' # note: this example will download 2 files of approx 10MB in total
#'
#' # store the downloaded files in the following directory. Here, the temporary file
#' # directory is used. Alternatively, consider storing this data in a more permanent location.
#' # e.g. output_dir="~/data/go" on unix systems or output_dir="C:/data/go" on Windows
#' output_dir = tempdir()
#'
#' # download data files with GO annotations (note that the release/date is printed to console)
#' # these are then parsed with the load_genesets_go_fromfile() function
#' # if the files are already available at output_dir, these are used and download is skipped
#' genesets_asis = download_genesets_goatrepo(output_dir)
#'
#' ### for a basic example on how to use the data obtain here,
#' ### refer to the example included at function documentation of: test_genesets()
#' }
#' @param output_dir full path to the directory where the downloaded files should be stored. Directory is created if it does not exist.
#' e.g. `output_dir="~/data"` on unix systems, `output_dir="C:/data"` on Windows, or set to `output_dir=getwd()` to write output to the current working directory
#' @param type the type of genesets to download. Currently, only "GO" is supported (default)
#' @param version the dataset version. This must be a date in format YYYY-MM-DD (for example; "2024-01-01") OR be left empty (NA or empty string, the default) to automatically download the latest version.
#' @param ignore_cache boolean, set to TRUE to force re-download and ignore cached data, if any. Default: FALSE
#' @return result from respective geneset parser function. e.g. if parameter `type` was set to"GO" (default), this function returns the result of `load_genesets_go_fromfile()`. These data returned by this function is typically used as input for `filter_genesets()`, c.f. full example at documentation for test_genesets()
#' @export
download_genesets_goatrepo = function(output_dir, type = "GO", version = "", ignore_cache = FALSE) {
  stopifnot("parameter type must be a single string. The only supported option for now is 'GO' (default)" = length(type) == 1 && is.character(type) && !is.na(type) && type %in% c("GO"))
  stopifnot("parameter version must be a single string that represents a date (indicates a version/release), or be left empty in order to automatically download the latest version" =
              length(version) == 1 && (!is.na(version) || (is.character(version) && (version == "" || grepl("^\\d\\d\\d\\d\\-\\d\\d-\\d\\d$", version)))) )
  stopifnot("parameter output_dir must be a single string and represent a directory on your computer" = length(output_dir) == 1 && is.character(output_dir) && !is.na(output_dir))
  stopifnot("parameter ignore_cache must be a single boolean value" = length(ignore_cache) == 1 && ignore_cache %in% c(TRUE, FALSE))

  # create dir if it does not exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if(!dir.exists(output_dir)) {
      stop(paste0("Could not create the requested output directory: ", output_dir, "\nTry to provide an existing directory as parameter for 'output_dir'"))
    }
  }

  if(type == "GO") {
    # if no version is provided, first download an overview of all available versions and identify the latest release
    if(is.na(version) || version == "") {
      all_versions = available_genesets_goatrepo()
      version = sort(unique(all_versions$date), decreasing = TRUE)[1] # sort all dates, take largest value
      message(paste("the most recent available is version:", version))
    }

    file_g2g = sprintf("%s/gene2go_%s.gz", output_dir, version)
    file_obo = sprintf("%s/go_%s.obo.gz", output_dir, version)
    url_g2g = sprintf("https://github.com/ftwkoopmans/goat/raw/data/go/%s/gene2go_human_%s.gz", version, version)
    url_obo = sprintf("https://github.com/ftwkoopmans/goat/raw/data/go/%s/go_%s.obo.gz", version, version)
    any_download = FALSE
    sprintf_template_downloadfail = "failed to download %s and store it at %s\nMost likely causes are 1) the requested file/version does not exist (try the default parameter!) and 2) Internet connection issues (try to download the here mentioned URL by copy/pasting in your browser).\nPlease refer to the download_genesets_goatrepo() function documentation to learn how you can find available versions (besides the default parameter)"

    # attempt to download if not available on disk
    if(ignore_cache || !file.exists(file_g2g)) {
      message(paste("downloading", url_g2g, "..."))
      utils::download.file(url_g2g, file_g2g, mode = "wb")
      if(!file.exists(file_g2g)) {
        stop(sprintf(sprintf_template_downloadfail, url_g2g, file_g2g))
      }
      any_download = TRUE
    }

    # attempt to download if not available on disk
    if(ignore_cache || !file.exists(file_obo)) {
      message(paste("downloading", url_obo, "..."))
      utils::download.file(url_obo, file_obo, mode = "wb")
      if(!file.exists(file_obo)) {
        stop(sprintf(sprintf_template_downloadfail, url_obo, file_obo))
      }
      any_download = TRUE
    }

    if(any_download) {
      message(paste("downloaded geneset files were stored at:", output_dir))
    } else {
      message(paste("cached geneset files were retrieved from:", output_dir))
    }

    return(load_genesets_go_fromfile(file_gene2go = file_g2g, file_goobo = file_obo))
  }

  # ... other types may be added in the future

}



#' Discover available geneset collections from the GOAT GitHub repository
#'
#' @description
#' Yields a table that describes available gene sets. You can download the actual data for the described versions with `download_genesets_goatrepo()`.
#' @export
available_genesets_goatrepo = function() {
  url = "https://github.com/ftwkoopmans/goat/raw/data/go/index.tsv"
  message(paste("downloading", url, "..."))
  x = data.table::fread(url, colClasses = "character", data.table = FALSE)
  if(length(x) < 2 || !is.data.frame(x) || nrow(x) == 0) {
    stop(paste("failed to download", url))
  }
  if(ncol(x) != 2 || !all(colnames(x) %in% c("date", "file"))) {
    print(x)
    stop("downloaded table should contain exactly 2 columns, with column names 'date' and 'file'")
  }
  return(tibble::as_tibble(x))
}
