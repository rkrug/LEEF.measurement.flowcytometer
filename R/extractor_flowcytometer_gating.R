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

  gates_coordinates <- utils::read.csv(file.path(input, "flowcytometer", "gates_coordinates.csv"))

  # function to gate each plate ---------------------------------------------

  gate_plate <- function(
    plate,
    input,
    output
  )  {


    #  read data files ----------------------------------------------------------


    flow.data <- readRDS(
      file = file.path(file.path(output, "flowcytometer"), paste0("flowcytometer_ungated.", plate, ".rds"))
    )
    fsa <- readRDS(
      file = file.path(file.path(output, "flowcytometer"), paste0("flowcytometer_fsa_ungated.", plate, ".rds"))
    )

    #############################################################
    # <<<< END SCRIPT   #########################################
    #############################################################

    # Apply the gating --------------------------------------------------------

    #RL: defining gates

    gates <- list()

    # bacteria gate
    polyGate_bacteria <- as.matrix(gates_coordinates[1:4, 1:2])
    colnames(polyGate_bacteria) <- c("FL1-A", "FL3-A")
    bacteria_gate <- flowCore::polygonGate(filterId = "Bacteria", .gate = polyGate_bacteria)

    # gate for different size classes of bacteria
    LNA_coordinates <- as.matrix(gates_coordinates[, 3])
    LNA_coordinates <- na.omit(LNA_coordinates)
    colnames(LNA_coordinates) <- c("FL1-A")

    MNA_coordinates <- as.matrix(gates_coordinates[, 4])
    MNA_coordinates <- na.omit(MNA_coordinates)
    colnames(MNA_coordinates) <- c("FL1-A")

    HNA_coordinates <- as.matrix(gates_coordinates[, 5])
    HNA_coordinates <- na.omit(HNA_coordinates)
    colnames(HNA_coordinates) <- c("FL1-A")

    rg_LNA <- flowCore::rectangleGate("FL1-A" = LNA_coordinates, filterId = "LNA")
    rg_MNA <- flowCore::rectangleGate("FL1-A" = MNA_coordinates, filterId = "MNA")
    rg_HNA <- flowCore::rectangleGate("FL1-A" = HNA_coordinates, filterId = "HNA")

    gates$bacteria <- list(
      bacteria_gate = bacteria_gate,
      rg_LNA = rg_LNA,
      rg_MNA = rg_MNA,
      rg_HNA = rg_HNA
    )

    # algae gate
    polyGate_algae <- as.matrix(gates_coordinates[1:4, 6:7])
    colnames(polyGate_algae) <- c("FL1-A", "FL4-A")
    algae_gate <- flowCore::polygonGate(filterId = "Algae", .gate = polyGate_algae)

    gates$algae <- list(
      algae_gate = algae_gate
    )

    # #----- HAVING A LOOK AT THE GATING -----#

    # i <- 1
    # flowViz::xyplot(`FL3-A` ~ `FL1-A`, data = fsa[[i]], filter = bacteria_gate)
    # # fsa[[i]]
    # flowViz::densityplot(~ `FL1-A`, data = fsa[[i]], filter = rg_LNA)
    # # fsa[[i]]
    # flowViz::xyplot(`FL4-A` ~ `FL1-A`, data = fsa[[i]], filter = algae_gate)
    # # fsa[[i]]

    # ABUNDANCE DYNAMICS ------------------------------------------------------

    # applying filter to whole flowSet

    result <- flowCore::filter(fsa, bacteria_gate)

    # extract absolute counts
    l <- lapply(result, flowCore::summary)

    counts <- sapply(
      l,
      function(i) {
        i$true
      }
    )


    flow.data$bacteria_counts <- counts
    flow.data$bacteria_density_perml <- flow.data$bacteria_counts * 1000000 / flow.data$volume * flow.data$dilution_factor
    flow.data$sample_letter <- substr(x = flow.data$sample, start = 1, stop = 1)
    flow.data$sample_number <- as.numeric(substr(x = flow.data$sample, start = 2, stop = 3))
    flow.data$date <- format(as.Date(flow.data$date, "%d-%b-%Y"), "%Y-%m-%d")
    flow.data <- flow.data[order(flow.data$date, flow.data$sample_letter, flow.data$sample_number), ]

    # subset data based on gate for bacteria
    subset.bacteria <- flowCore::Subset(fsa, bacteria_gate)

    # applying filter to bacteria to get the three bacteria populations
    LNA <- flowCore::filter(subset.bacteria, rg_LNA)
    MNA <- flowCore::filter(subset.bacteria, rg_MNA)
    HNA <- flowCore::filter(subset.bacteria, rg_HNA)

    # extract absolute counts
    l_LNA <- lapply(LNA, flowCore::summary)
    l_MNA <- lapply(MNA, flowCore::summary)
    l_HNA <- lapply(HNA, flowCore::summary)

    counts_LNA <- sapply(
      l_LNA,
      function(i) {
        i$true
      }
    )
    counts_MNA <- sapply(
      l_MNA,
      function(i) {
        i$true
      }
    )
    counts_HNA <- sapply(
      l_HNA,
      function(i) {
        i$true
      }
    )

    flow.data$LNA_counts <- counts_LNA
    flow.data$MNA_counts <- counts_MNA
    flow.data$HNA_counts <- counts_HNA

    flow.data$LNA_perml <- flow.data$LNA_counts * 1000000 / flow.data$volume * flow.data$dilution_factor
    flow.data$MNA_perml <- flow.data$MNA_counts * 1000000 / flow.data$volume * flow.data$dilution_factor
    flow.data$HNA_perml <- flow.data$HNA_counts * 1000000 / flow.data$volume * flow.data$dilution_factor

    # get the algae
    algae <- flowCore::filter(fsa, algae_gate)

    # extract absolute counts
    l_algae <- lapply(
      algae,
      flowCore::summary
    )
    counts_algae <- sapply(
      l_algae,
      function(i) {
        i$true
      }
    )

    flow.data$algae_counts <- counts_algae
    flow.data$algae_perml <- flow.data$algae_counts * 1000000 / flow.data$volume * flow.data$dilution_factor

    ### Change to long format ###
    flow.data1 <- flow.data %>%
      dplyr::select(-LNA_perml, -MNA_perml, -HNA_perml, -algae_perml, -bacteria_density_perml) %>%
      dplyr::rename(LNA = LNA_counts, MNA = MNA_counts, HNA = HNA_counts,
                    algae = algae_counts, bacteria = bacteria_counts) %>%
      tidyr::pivot_longer(names_to = "name", values_to = "count",
                          cols = c(LNA, MNA, HNA, algae, bacteria),
                          names_transform = list())

    flow.data2 <- flow.data %>%
      dplyr::select(-LNA_counts, -MNA_counts, -HNA_counts, -algae_counts, -bacteria_counts) %>%
      dplyr::rename(LNA = LNA_perml, MNA = MNA_perml, HNA = HNA_perml,
                    algae = algae_perml, bacteria = bacteria_density_perml) %>%
      tidyr::pivot_longer(names_to = "name", values_to = "density",
                          cols = c(LNA, MNA, HNA, algae, bacteria),
                          names_transform = list())

    flow.data <- plyr::join(flow.data1, flow.data2)

    #############################################################
    # >>>> END SCRIPT   #########################################
    #############################################################

    # Rename 'name' column to 'species' for consistency -----------------------

    names(flow.data)[which(names(flow.data) == "name")] <- "species"

    # SAVE --------------------------------------------------------------------

    saveRDS(
      flow.data,
      file = file.path(add_path, paste0("flowcytometer_density.", plate, ".rds"))
    )

    saveRDS(
      gates,
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
    function(x){
      plate <-strsplit(x, split = "\\.")[[1]][[2]]
      return(plate)
    }
  )
  plates <- unique(plates)

  lapply(
    plates,
    gate_plate,
    input = input,
    output = output
  )

  # Finalise ----------------------------------------------------------------

  flow.data <- NULL
  for (plate in plates){
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

# Finalize ----------------------------------------------------------------

  unlink(processing)
  message("   done")
  message("########################################################")

  invisible(TRUE)
}
