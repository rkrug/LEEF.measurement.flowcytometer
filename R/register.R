#' Register the processing of flowcytometer data in the LEEF package
#'
#' @return invisibly \code{TRUE} when completed successful
#'

#' @export
#'
register <- function() {
  if (is.null(system.file(package = "LEEF.2"))) {
    stop("This function requres the package to be installed!")
  }

  LEEF.2::add_pre_processor( pre_processor_flowcytometer )
  LEEF.2::add_extractor( extractor_flowcytometer)

  invisible(TRUE)
}

