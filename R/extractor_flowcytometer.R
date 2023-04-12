#' Extractor flowcytometer data
#'
#'
#' This function is extracting data to be added to the database (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#' @param min_FSC.A numeric. If \code{!NULL}, \code{FSA.A <= min_FSC.A} will be fitered out by using
#'   a rectangular filter
#'   \code{flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))}
#' @param use_H if \code{TRUE}, gating will be done using \code{height}, otherwie \code{area}
#'
#' @return invisibly \code{TRUE} when completed successful
#'
#' @export
#'
extractor_flowcytometer <- function(
  input,
  output,
  raw = FALSE,
  min_FSC.A = NULL,
  use_H = FALSE
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
    }
  )
  on.exit({
      if (file.exists(processing)) {
        unlink(processing)
        file.create(error)
      }
    }
  )

  ##

  extractor_flowcytometer_preparation(
    input,
    output,
    raw = raw,
    min_FSC.A = min_FSC.A
  )
  extractor_flowcytometer_density(
    input,
    output,
    min_FSC.A = min_FSC.A,
    use_H = use_H
  )

# Finalize ----------------------------------------------------------------

  unlink(processing)
  message("done")
  message("########################################################")

  invisible(TRUE)
}
