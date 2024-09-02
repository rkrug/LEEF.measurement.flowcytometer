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
extractor_flowcytometer_density <- function(
    input = NA,
    output = ".",
    min_FSC.A = NULL,
    use_H = FALSE,
    gates_coordinates = NULL,
    fsa = NULL,
    flow.data = NULL,
    dens_back = FALSE) {
  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  log.file <- tempfile(pattern="flowcytometer.", tmpdir = add_path, fileext = ".log")
  loggit::set_logfile(log.file)

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
  })
  on.exit({
    if (file.exists(processing)) {
      unlink(processing)
      file.create(error)
    }
  })

  ##

  #############################################################
  #############################################################




  #############################################################
  #############################################################

  # function to gate each plate ---------------------------------------------

  gate_plate <- function(
      plate,
      input,
      output,
      use_H,
      gates_coordinates,
      fsa = NULL,
      flow.data = NULL) {
    if (is.null(gates_coordinates)) {
      gates_coordinates <- utils::read.csv(file.path(input, "flowcytometer", "gates_coordinates.csv"))
    }

    if (is.null(fsa)) {
      fsa <- readRDS(file.path(output, "flowcytometer", paste0("flowcytometer_fsa_ungated.", plate, ".rds")))
    }

    if (is.null(flow.data)) {
      flow.data <- readRDS(file.path(output, "flowcytometer", paste0("flowcytometer_ungated.", plate, ".rds")))
    }

    if (use_H) {
      gated <- dens_H(
        gates_coordinates = gates_coordinates,
        fsa = fsa,
        flow.data = flow.data
      )
    } else {
      gated <- dens(
        gates_coordinates = gates_coordinates,
        fsa = fsa,
        flow.data = flow.data
      )
    }
    # SAVE --------------------------------------------------------------------

    saveRDS(
      gated$flow.data,
      file = file.path(add_path, paste0("flowcytometer_density.", plate, ".rds"))
    )

    saveRDS(
      gated$gates,
      file = file.path(add_path, paste0("flowcytometer_gates.", plate, ".rds"))
    )
  }


  # Do the gating for all plates --------------------------------------------


  plates <- grep(
    "p_",
    list.files(
      file.path(output, "flowcytometer"),
      full.names = FALSE
    ),
    value = TRUE
  )
  plates <- sapply(
    plates,
    function(x) {
      plate <- strsplit(x, split = "\\.")[[1]][[2]]
      return(plate)
    }
  )
  plates <- unique(plates)

  lapply(
    plates,
    gate_plate,
    input = input,
    output = output,
    use_H = use_H,
    gates_coordinates,
    fsa = fsa,
    flow.data = flow.data
  )

  # Finalise ----------------------------------------------------------------

  if (!dens_back) {
    flow.data <- NULL
    for (plate in plates) {
      fdp <- readRDS(file.path(add_path, paste0("flowcytometer_density.", plate, ".rds")))
      flow.data <- rbind(flow.data, fdp)
    }
    utils::write.csv(
      flow.data,
      file = file.path(add_path, paste0("flowcytometer_density.csv")),
      row.names = FALSE
    )
    unlink(list.files(add_path, "flowcytometer_density\\.p_.\\.rds", full.names = TRUE))

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
    density <- NULL
    for (plate in plates) {
      fdp <- readRDS(file.path(add_path, paste0("flowcytometer_density.", plate, ".rds")))
      density <- rbind(density, fdp)
    }
    return(density)
  } else {
    invisible(TRUE)
  }
}
