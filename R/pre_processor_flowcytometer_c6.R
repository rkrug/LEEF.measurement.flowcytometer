#' Preprocessor flowcytometer data FOR c6 FILES!!!
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
#'
#' @export
#'
pre_processor_flowcytometer_c6 <- function(
  input,
  output
) {
  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  loggit::set_logfile(file.path(add_path, "flowcytometer.log"))

  ##
  message("########################################################")
  message("Processing flowcytometer...")
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

  oldwd <- getwd()
  setwd( file.path( tmp ) )

  cmd <- "python"
  arguments <- system.file(package = "LEEF.measurement.flowcytometer", "accuri2fcs", "accuri2fcs", "accuri2fcs.py" )
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

  setwd(oldwd)

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

  message("done")
  message("########################################################")

  invisible(TRUE)
}
