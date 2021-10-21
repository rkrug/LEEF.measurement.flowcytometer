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

  if ( length( list.files( file.path(input, "flowcytometer") ) ) == 0 ) {
    message("Empty or missing flowcytometer directory - nothing to do.")
    message("done")
    message("########################################################")
    return(invisible(TRUE))
  }

  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  loggit::set_logfile(file.path(add_path, "flowcytometer.log"))

  ##
  message("########################################################")
  message("Processing flowcytometer...")


  dir.create(
    file.path(output, "flowcytometer"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  ##


  file.copy(
    file.path( input, "..", "00.general.parameter", "." ),
    file.path( output, "flowcytometer" ),
    recursive = TRUE,
    overwrite = TRUE
  )

  file.copy(
    from = file.path(input, "flowcytometer", "."),
    to = file.path(output, "flowcytometer"),
    recursive = TRUE
  )

  ##


  message("done")
  message("########################################################")
  ##
  invisible(TRUE)
}
