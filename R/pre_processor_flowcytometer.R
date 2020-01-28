#' Preprocessor flowcytometer data
#'
#' Convert all \code{.c6} files in \code{flowcytometrie} folder to \code{.fcs} files and delete \code{.c6} file.
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
pre_processor_flowcytometer <- function( input, output ) {
  ##
  oldwd <- getwd()
  on.exit(
    setwd(oldwd)
  )
  ##

  message("\n########################################################\n")
  message("\nProcessing flowcytometer...\n")
  setwd( file.path( input, "flowcytometer" ) )
  cmd <- "python"
  arguments <- file.path( system.file(package = utils::packageName()), "tools", "accuri2fcs", "accuri2fcs", "accuri2fcs.py" )
  system2(
    command = cmd,
    args = arguments
  )
  unlink("*.c6")
  fcs <- list.files(
    path = file.path( input, "flowcytometer", "fcs" ),
    full.names = TRUE
  )
  file.rename(
    from = fcs,
    to = gsub( "/fcs/", "/", fcs )
  )
  unlink( file.path( input, "flowcytometer", "fcs" ), recursive = TRUE )
  message("done\n")
  message("\n########################################################\n")

  invisible(TRUE)
}
