
#' ASCII logo for this package
#'
#' @return package logo as a string
#' @export
goat_logo = function() {
  '          )_)
  (\\_____\'\\_\\
  |      |  "
  |""""""|
'
}



#' Return goat package version as a string
#'
#' simple wrapper around utils::packageVersion()
#' @return package version as a string
#' @export
goat_version = function() {
  as.character(utils::packageVersion("goat"))
}



#' Print package version and logo to console
#'
#' @return prints to console without returning a value
#' @export
goat_print_version = function() {
  message(sprintf("%s\nGOAT version %s", goat_logo(), goat_version()))
}



#' return a prettyprint string of all length 1 parameters that are string/numeric/logical/NA
#'
#' @examples
#' parameters_prettyprint_length1(
#'   test1=1:2, test2=matrix(1:4,2,2), test3=data.frame(a=1),
#'   test4=c(a=1), test5=1, test6="a", test7=NA, test8=Inf
#' )
#' @param ... arbitrary set of parameters
#' @noRd
parameters_prettyprint_length1 = function(...) {
  arguments = list(...)
  if(length(arguments) == 0 || !is.list(arguments)) {
    return()
  }

  result = NULL
  for(n in names(arguments)) {
    x = arguments[[n]]
    if(length(x) == 1) {
      if(is.character(x) && !is.na(x)) {
        result = c(result, paste0(n,"='", x, "'"))
      } else if(is.na(x) || is.numeric(x) || is.logical(x)) {
        result = c(result, paste0(n,"=", x))
      }
    }
  }

  paste(result, collapse = ", ") # when result is NULL, return empty string
}



#' throw error if R package is unavailable
#'
#' @param pkg R package name
#' @param msg function name / reference for user
#' @noRd
check_dependency = function(pkg, msg) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0(
      "An optional dependency for ",
      msg,
      " is not installed; R package '",
      pkg,
      "' is not available. For convenience, you may use the following command to install all dependencies (including optional) for the 'goat' R package; pak::pkg_install('ftwkoopmans/goat', dependencies = TRUE)"
    ), call. = FALSE)
  }
}



#' generate colours analogous to ggplot's default palette
#'
#' https://stackoverflow.com/a/8197703
#'
#' @param n number of colors
#' @return a color code (string)
#' @export
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}



#' naively lighten a color by mixing in white
#'
#' @param color input colors
#' @param frac fraction of white; >0 and <1
#' @return adjusted value for input `color`
#' @export
lighten_color = function(color, frac = 0.1){
  stopifnot(frac > 0 & frac < 1)
  sapply(color, function(x) grDevices::colorRampPalette(c(x, "white"))(100)[ceiling(frac * 100)], simplify = TRUE, USE.NAMES = FALSE)
}



#' naively darken a color by mixing in black
#'
#' @param color input colors
#' @param frac fraction of black; >0 and <1
#' @return adjusted value for input `color`
#' @export
darken_color = function(color, frac = 0.1){
  stopifnot(frac > 0 & frac < 1)
  sapply(color, function(x) grDevices::colorRampPalette(c(x, "black"))(100)[ceiling(frac * 100)], simplify = TRUE, USE.NAMES = FALSE)
}



#' check if a string is a valid R color-code
#'
#' @param x input string
#' @noRd
isvalid_color = function(x) {
  length(x) == 1 && !is.na(x) && is.character(x) && suppressWarnings(tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE))
}



#' simple string truncation
#'
#' replacement for stringr::trunc() so we don't need a package dependency for just 1 function (our code was adapter therefrom)
#' @param string string that should be truncated
#' @param width desired max length
#' @param trim_left instead of right trunc (default), do left instead
#' @return truncated variant of input `string`
#' @export
string_trunc_right = function(string, width, trim_left = FALSE) {
  N = nchar(string)
  too_long = !is.na(string) & N > 3 & N > width # hardcoded min string length @ 4
  if(any(too_long)) {
    if(trim_left) {
      string[too_long] = paste0("...", substr(string[too_long], N - width - 3, N))
    } else {
      string[too_long] = paste0(substr(string[too_long], 1, width - 3), "...")
    }
  }
  string
}
