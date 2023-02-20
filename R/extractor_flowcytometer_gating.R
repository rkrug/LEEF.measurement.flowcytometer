#' Extractor flowcytometer data
#'
#'
#' This function is extracting data to be added to the database
#' (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#'
#' @return invisibly \code{TRUE} when completed successful
#'
#' @importFrom flowCore read.flowSet pData phenoData exprs logTransform truncateTransform transform
#' @importFrom flowCore rectangleGate polygonGate Subset
#' @importFrom yaml read_yaml
#' @importFrom stats setNames
#' @importFrom utils read.csv write.csv
#' @importFrom plyr join
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @import loggit
#' @export
#'
extractor_flowcytometer_gating <- function(
    input,
    output
) {
  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  loggit::set_logfile(file.path(add_path, "flowcytometer.log"))

  message("########################################################")
  message("   gating flowcytometer...")

  ##
  suppressWarnings({
    processing <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("EXTRACTING.FLOWCYTOMETER.GATING", ".PROCESSING")
    )
    error <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("ERROR.EXTRACTING.FLOWCYTOMETER.GATING", ".ERROR")
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



  #############################################################
  #############################################################


  flow.data <- gating(
    gates_coordinates = utils::read.csv(file.path(input, "flowcytometer", "gates_coordinates.csv")),
    fsa = readRDS(file.path(output, "flowcytometer", "flowcytometer_fsa_ungated.rds")),
    flow.data = utils::read.csv(file.path(output, "flowcytometer", "flowcytometer_ungated.csv"))
  )$flow.data



  # SAVE --------------------------------------------------------------------

  utils::write.csv(
    flow.data,
    file = file.path(add_path, "flowcytometer_density.csv"),
    row.names = FALSE
  )
  to_copy <- grep(
    list.files(
      file.path(input, "flowcytometer"),
      full.names = TRUE
    ),
    pattern = "\\.ciplus$",
    invert = TRUE,
    value = TRUE
  )
  file.copy(
    from = to_copy,
    to = file.path(output, "flowcytometer", "")
  )

  # Finalize ----------------------------------------------------------------

  unlink(processing)
  message("   done")
  message("########################################################")

  invisible(TRUE)
}
