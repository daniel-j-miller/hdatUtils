#' Title
#'
#' @param dat a dataframe imported with col name attrs
#'
#' @returns vector of column names and attrs
#'
#' @importFrom utils writeClipboard
#' @export
get_dta_cols <- function(dat) {

  tmp_cols <- vapply(names(dat),
                     function(x) paste(x,"~",attr(dat[[x]], "label")),character(1))
  writeClipboard(tmp_cols)
  tmp_cols
}
