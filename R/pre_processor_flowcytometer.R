#' Preprocessor flowcytometer data
#'
#' Just copy all files
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#'
#' @return invisibly \code{TRUE} when completed successful
#'
#' @export
#'
pre_processor_flowcytometer <- function(
  input,
  output
) {
  ##
  message("\n########################################################\n")
  message("\nProcessing flowcytometer...\n")
  ##
  dir.create(
    file.path(output, "flowcytometer"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  file.copy(
    from = file.path(input, "flowcytometer", "."),
    to = file.path(output, "flowcytometer"),
    recursive = TRUE
  )
  file.copy(
    from = file.path(input, "sample_metadata.yml"),
    to = file.path(output, "sample_metadata.yml")
  )

  ##
  message("done\n")
  message("\n########################################################\n")
  ##
  invisible(TRUE)
}
