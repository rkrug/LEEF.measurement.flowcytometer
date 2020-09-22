#' Check if data in input folder is OK and move to raw data folder
#'
#' @param input The folder, where a folder \code{flowcytometer} is located which
#'   contains the new files.
#' @param output A folder, which contains a subfolder called \code{flowcytometer}, i.e.
#'   the usually the raw data folder, into which the fioles will be moved to.
#'
#' @return a \code{list} which contains the individual reseults for each file.
#'   \code{TRUE} if moved, \code{FALSE} if an error occured. Details of the eror
#'   re in the error files in the \code{input/flowcytometer} directory.
#' @importFrom parallel mclapply
#' @export
#'
add_new_data <- function(input, output) {
  ##
  dir.create(
    file.path(output, "flowcytometer"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  # Copy ALL other files ----------------------------------------------------

#   others <- grep(
#     list.files(
#       path = input,
#       full.names = TRUE
#     ),
#     pattern='.cxd',
#     invert=TRUE,
#     value=TRUE
#   )
#   file.copy(
#     from = others,
#     to = file.path(output, "flowcytometer"),
#     overwrite = TRUE
#   )
#   unlink( others )

  # Check and move folder ------------------------------------------------------

  files <- list.files(
    path = input,
    pattern = "\\.c6$",
    full.names = FALSE
  )

  ##
  ok <- parallel::mclapply(
    files,
    function(f) {
      processing <- file.path(input, paste0("CHECKING.", f, ".CHECKING"))
      error <- file.path(input, paste0("ERROR.", f, ".txt"))

      on.exit(
        {
          if (file.exists(processing)) {
            unlink(processing)
            capture.output(print(result), file = error)
          }
        }
      )
      ##
      file.create( processing )
      ##
      message("checking ", f)
      result <- list(
        ok = TRUE
      )

      # Check if file exist ----------------------------------------------------------


      result$exists <- file.exists(file.path(input, f))

      result$ok <- all(unlist(result))

      if ( result$ok ) {
        file.copy(
          from = file.path(input, f),
          to = file.path(output, "flowcytometer"),
          recursive = FALSE,
          overwrite = TRUE
        )
        unlink( file.path(input, f), recursive = TRUE )
        unlink(processing)
      }
      return(result)
    }
  )
  names(ok) <- files
  return(ok)
}
