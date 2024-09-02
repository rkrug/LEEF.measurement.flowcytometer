#' Extractor flowcytometer data
#'
#'
#' This function is extracting data to be added to the database (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#'
#' @return invisibly \code{TRUE} when completed successful
#'
#' @export
#'
extractor_flowcytometer <- function(
    input,
    output) {
  if (length(list.files(file.path(input, "flowcytometer"))) == 0) {
    message("Empty or missing flowcytometer directory - nothing to do.")
    message("done")
    message("########################################################")
    return(invisible(TRUE))
  }

  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  log.file <- tempfile(pattern="flowcytometer.", tmpdir = add_path, fileext = ".log")
  loggit::set_logfile(log.file)

  message("########################################################")
  message("Extracting flowcytometer...")

  ##
  suppressWarnings({
    processing <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("EXTRACTING.FLOWCYTOMETER", ".PROCESSING")
    )
    error <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("ERROR.EXTRACTING.FLOWCYTOMETER", ".ERROR")
    )
    file.create(processing)
  })
  on.exit({
    if (file.exists(processing)) {
      unlink(processing)
      file.create(error)
    }
  })

  ##

  extractor_flowcytometer_preparation(input, output)
  extractor_flowcytometer_density(input, output)
  extractor_flowcytometer_traits(input, output)
  ## Traits

  # Finalize ----------------------------------------------------------------

  unlink(processing)
  message("done")
  message("########################################################")

  invisible(TRUE)
}
