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
#' @importFrom flowCore read.flowSet pData phenoData exprs logTransform truncateTransform transform rectangleGate Subset
#' @importFrom  yaml read_yaml
#' @importFrom  stats setNames
#' @importFrom utils write.csv
#' @export
#'
extractor_flowcytometer <- function(
  input,
  output
) {
  message("\n########################################################\n")
  message("Extracting flowcytometer gating...\n")

  ##
  suppressWarnings(
    {
      processing <- file.path(normalizePath(output), "flowcytometer", paste0("EXTRACTING.FLOWCYTOMETER.GATING", ".PROCESSING"))
      error <- file.path(normalizePath(output), "flowcytometer", paste0("ERROR.EXTRACTING.FLOWCYTOMETER.GATING", ".ERROR"))
      file.create( processing )
    }
  )
  on.exit(
    {
      if (file.exists(processing)) {
        unlink(processing)
        file.create(error)
      }
    }
  )

  ##

  gates_coordinates <- read.csv(file.path( input, "flowcytometer", "gates_coordinates.csv" ))

  fs <- file.path( output, "flowcytometer", "flowcytometer_ungated.csv" )


  # Based on flowcyt_1_c6_to_RData.R ----------------------------------------
  # Converting the Flowcytometer Output of bacterial abundances into a usable data frame
  # David Inauen, 19.06.2017


  # Apply the gating --------------------------------------------------------

  #RL: defining gates

  # bacteria gate
  polyGate_bacteria <- as.matrix(gates_coordinates[1:4, 1:2])
  colnames(polyGate_bacteria) <- c("FL1-A","FL3-A")
  bacteria_gate <- flowCore::polygonGate(filterId="Bacteria", .gate = polyGate_bacteria)

  # gate for different size classes of bacteria
  LNA_coordinates <- as.matrix(gates_coordinates[,3])
  LNA_coordinates <- na.omit(LNA_coordinates)
  colnames(LNA_coordinates) <- c("FL1-A")

  MNA_coordinates <- as.matrix(gates_coordinates[,4])
  MNA_coordinates <- na.omit(MNA_coordinates)
  colnames(MNA_coordinates) <- c("FL1-A")

  HNA_coordinates <- as.matrix(gates_coordinates[,5])
  HNA_coordinates <- na.omit(HNA_coordinates)
  colnames(HNA_coordinates) <- c("FL1-A")

  rg_LNA <- rectangleGate("FL1-A"= LNA_coordinates, filterId="LNA")
  rg_MNA <- rectangleGate("FL1-A"= MNA_coordinates, filterId="MNA")
  rg_HNA <- rectangleGate("FL1-A"= HNA_coordinates, filterId="HNA")


  # algae gate
  polyGate_algae <- as.matrix(gates_coordinates[1:4, 6:7])
  colnames(polyGate_algae) <- c("FL1-A","FL4-A")
  algae_gate <- flowCore::polygonGate(filterId="Algae", .gate= polyGate_algae)

# ABUNDANCE DYNAMICS ------------------------------------------------------

# Apply bacteria gate -----------------------------------------------------

   # applying filter to whole flowSet
   result <- flowCore::filter(fsa, bacteria_gate)

   # extract absolute counts
   l <- lapply(result, flowCore::summary)
   # counts <- plyr::ldply(
   #   lapply(
   #     l,
   #     function(i) {
   #       i$true
   #     }
   #   )
   # )

   counts <- sapply(
     l,
     function(i) {
       i$true
     }
   )
   flow.data[["bacteria_counts"]] <- counts
   flow.data[["bacteria_density_perml"]] <- flow.data[["bacteria_counts"]] * 1000000/flow.data[["volume"]] * flow.data[["dilution_factor"]]
   flow.data$sample_letter <- substr(
     x = flow.data$sample,
     start = 1,
     stop = 1
   )
  flow.data[["sample_number"]] <- as.numeric(
    substr(
      x = flow.data[["sample"]],
      start = 2,
      stop = 3
    )
  )

  # this does not run
  # flow.data[["date"]] <- format(
  #   as.Date(
  #     flow.data$date,
  #     "%d-%b-%Y"
  #   ),
  #   "%Y-%m-%d"
  # )
  flow.data <- flow.data[ order(
#     flow.data$date,
    flow.data$sample_letter,
    flow.data$sample_number
  ),]


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


  flow.data[["LNA_counts"]] <- counts_LNA
  flow.data[["MNA_counts"]] <- counts_MNA
  flow.data[["HNA_counts"]] <- counts_HNA

  flow.data[["LNA_perml"]] <- flow.data[["LNA_counts"]] * 1000000/flow.data[["volume"]] * flow.data[["dilution_factor"]]
  flow.data[["MNA_perml"]] <- flow.data[["MNA_counts"]] * 1000000/flow.data[["volume"]] * flow.data[["dilution_factor"]]
  flow.data[["HNA_perml"]] <- flow.data[["HNA_counts"]] * 1000000/flow.data[["volume"]] * flow.data[["dilution_factor"]]


# Apply alagae gate -------------------------------------------------------

  # get the algae
  algae <- flowCore::filter(fsa, algae_gate)

  # extract absolute counts
  l_algae <- lapply(algae, flowCore::summary)

  counts_algae <- sapply(
    l_algae,
    function(i) {
      i$true
    }
  )

  flow.data[["algae_counts"]] <- counts_algae
  flow.data[["algae_perml"]] <- flow.data[["algae_counts"]] * 1000000/flow.data[["volume"]] * flow.data[["dilution_factor"]]

# SAVE --------------------------------------------------------------------

  add_path <- file.path( output, "flowcytometer" )
  dir.create( add_path, recursive = TRUE, showWarnings = FALSE )

  timestamp <- yaml::read_yaml(file.path(input,  "flowcytometer", "sample_metadata.yml"))$timestamp

  flow.data <- cbind(
  	timestamp = timestamp,
  	flow.data
  )

  utils::write.csv(
    flow.data,
    file = file.path(add_path, "flowcytometer.csv"),
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
  message("done\n")
  message("\n########################################################\n")

  invisible(TRUE)
}
