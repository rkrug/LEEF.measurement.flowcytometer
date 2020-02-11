#' Preprocessor flowcytometer data
#'
#' Convert all \code{.c6} files in \code{flowcytometrie} folder to \code{.fcs} files in output folder
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#'
#' @return invisibly \code{TRUE} when completed successful
#'
#' @importFrom R.utils bzip2
#' @importFrom tiff readTIFF writeTIFF
#' @importFrom parallel mclapply detectCores
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
  oldwd <- getwd()
  on.exit(
    setwd(oldwd)
  )
  ##
  tmp <- tempfile()
  dir.create(tmp)
  ##
  file.copy(
    from = file.path(input, "flowcytometer", "."),
    to = tmp,
    recursive = TRUE
  )


# Convert from .c6 to .fcs ------------------------------------------------


  setwd( file.path( tmp ) )
  cmd <- "python"
  arguments <- system.file(package = "LEEF.measurement.flowcytometer", "tools", "accuri2fcs", "accuri2fcs", "accuri2fcs.py" )
  system2(
    command = cmd,
    args = arguments,
    stdout = ifelse(
      options()$LEEF.measurement.flowcytometer$debug,
      "",
      FALSE
    )
  )
  unlink("*.c6")
  ##
  fcs <- list.files(
    path = file.path( tmp, "fcs" ),
    full.names = TRUE
  )
  file.rename(
    from = fcs,
    to = gsub( "/fcs/", "/", fcs )
  )
  unlink( file.path( tmp, "fcs" ), recursive = TRUE )


# Copy to final output ----------------------------------------------------


  dir.create(
    file.path(output, "flowcytometer"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  file.copy(
    from = file.path(tmp, "."),
    to = file.path(output, "flowcytometer"),
    recursive = TRUE
  )

  message("done\n")
  message("\n########################################################\n")

  invisible(TRUE)
}
