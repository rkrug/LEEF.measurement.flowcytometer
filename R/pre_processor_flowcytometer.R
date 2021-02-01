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
  processing <- file.path(normalizePath(output), "bemovi", paste0("PRE-PROCESSING.FLOWCYTOMETER", ".PROCESSING"))
  error <- file.path(normalizePath(output), "bemovi", paste0("ERROR.PRE-PROCESSING.FLOWCYTOMETER", ".ERROR"))
  on.exit(
    {
      if (file.exists(processing)) {
        unlink(processing)
        file.create(error)
      }
    }
  )
  file.create( processing )
  ##


  ##
  dir.create(
    file.path(output, "flowcytometer"),
    recursive = TRUE,
    showWarnings = FALSE
  )
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

  unlink(processing)
  message("done\n")
  message("\n########################################################\n")
  ##
  invisible(TRUE)
}
