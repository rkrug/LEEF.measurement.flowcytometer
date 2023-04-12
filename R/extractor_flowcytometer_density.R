#' Extractor flowcytometer density
#'
#'
#' This function is extracting data to be added to the database
#' (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param output directory to which to write the data
#' @param min_FSC.A numeric. If \code{!NULL}, \code{FSA.A <= min_FSC.A} will be fitered out by using
#'   a rectangular filter
#'   \code{flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))}
#' @param use_H if \code{TRUE}, gating will be done using \code{height}, otherwie \code{area}
#' @param gates_coordinates if \code{NULL}, \code{gates_coordinates} will be read in, otherwise the \code{gates_coordinates}
#' @param fsa if \code{NULL}, \code{fsa} will be read in, otherwise the \code{fsa}
#' @param flow.data if \code{NULL}, \code{flow.data} will be read in, otherwise the \code{flow.data}
#' @param dens_back if \code{FALSE}, density will be saved in "flowcytometer_density.csv", otherwise it will be returned
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
extractor_flowcytometer_density <- function(
    input = NULL,
    output = ".",
    min_FSC.A = NULL,
    use_H = FALSE,
    gates_coordinates = NULL,
    fsa = NULL,
    flow.data = NULL,
    dens_back = FALSE
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

  if (is.null(gates_coordinates)){
    gates_coordinates <- utils::read.csv(file.path(input, "flowcytometer", "gates_coordinates.csv"))
  }

  if (is.null(fsa)){
    fsa <- readRDS(file.path(output, "flowcytometer", "flowcytometer_fsa_ungated.rds"))
  }

  if (is.null(flow.data)){
    flow.data <-  utils::read.csv(file.path(output, "flowcytometer", "flowcytometer_ungated.csv"))
  }

  if (use_H) {
    flow.data <- dens_H(
      gates_coordinates = gates_coordinates,
      fsa = fsa,
      flow.data = flow.data,
      min_FSC.A = min_FSC.A
    )$flow.data
  } else {
    flow.data <- dens(
      gates_coordinates = gates_coordinates,
      fsa = fsa,
      flow.data = flow.data,
      min_FSC.A = min_FSC.A
    )$flow.data
  }


  # SAVE --------------------------------------------------------------------

  if (!dens_back) {
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
  }

  # Finalize ----------------------------------------------------------------

  unlink(processing)
  message("   done")
  message("########################################################")

  if (dens_back) {
    return(flow.data)
  } else {
    invisible(TRUE)
  }
}
