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
extractor_flowcytometer_traits <- function(
  input,
  output
) {
  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  loggit::set_logfile(file.path(add_path, "flowcytometer.log"))

  message("########################################################")
  message("   extracting traits flowcytometer...")

  ##
  suppressWarnings({
      processing <- file.path(
        normalizePath(output), "flowcytometer",
          paste0("EXTRACTING.FLOWCYTOMETER.TRAITS", ".PROCESSING")
        )
      error <- file.path(
        normalizePath(output), "flowcytometer",
        paste0("ERROR.EXTRACTING.FLOWCYTOMETER.TRAITS", ".ERROR")
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

  trait_plate <- function(
    plate,
    input,
    output
  )  {
    #  read data files ----------------------------------------------------------

    message("   reading data ", plate, " ...")
    # flow.data <- readRDS(
    #   file = file.path(file.path(output, "flowcytometer"), paste0("flowcytometer_ungated.", plate, ".rds"))
    # )
    fsa <- readRDS(
      file = file.path(file.path(output, "flowcytometer"), paste0("flowcytometer_fsa_ungated.", plate, ".rds"))
    )


    gates <- readRDS(file.path(output, "flowcytometer", paste0("flowcytometer_gates.", plate, ".rds")))


    # The Bacteria ------------------------------------------------------------

    message("   gating bacteria ", plate, " ...")


    bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)

    # LNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_LNA)
    # MNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_MNA)
    # HNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_HNA)


    # The Algae ---------------------------------------------------------------


    message("   gating algae ", plate, " ...")

    algae_pop <- Subset(fsa,  gates$algae$algae_gate)


    # extraction function -----------------------------------------------------


    extr_traits <- function(
      pop
    ){
      traits <- flowCore::fsApply(
        pop,
        function(p){
          result <- list()
          result$sample <- unlist(flowCore::keyword(p, "$WELLID"))
          result$plate <- plate
          x <- exprs(p)
          if (nrow(x) > 0) {
            result <- suppressWarnings(
              data.frame(
                result,
                x
              )
            )
          } else {
            result <- NULL
          }
          return(result)
        }
      )

      traits <- do.call(rbind, traits)
      return(traits)
    }


    # SAVE --------------------------------------------------------------------

    message("   extracting bacteria traits ", plate, " ...")
    saveRDS(
      extr_traits(bacteria_pop),
      file = file.path(output, "flowcytometer", paste0("flowcytometer_traits_bacteria.", plate, ".rds"))
    )
    message("   extracting LNA traits ", plate, " ...")
    # saveRDS(
    #   extr_traits(LNA_pop),
    #   file = file.path(output, "flowcytometer", paste0("flowcytometer_traits_lna.", plate, ".rds"))
    # )
    # message("   extracting MNA traits ", plate, " ...")
    # saveRDS(
    #   extr_traits(MNA_pop),
    #   file = file.path(output, "flowcytometer", paste0("flowcytometer_traits_mna.", plate, ".rds"))
    # )
    # message("   extracting HNA traits ", plate, " ...")
    # saveRDS(
    #   extr_traits(HNA_pop),
    #   file = file.path(output, "flowcytometer", paste0("flowcytometer_traits_hna.", plate, ".rds"))
    # )
    message("   extracting algae traits ", plate, " ...")
    saveRDS(
      extr_traits(algae_pop),
      file = file.path(output, "flowcytometer", paste0("flowcytometer_traits_algae.", plate, ".rds"))
    )

  }


  # Do the trait extraction for all plates --------------------------------------------


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
    trait_plate,
    input = input,
    output = output
  )


  traits <- NULL
  #
  for (plate in plates){
    traits <- rbind(
      traits,
      readRDS(file.path(add_path, paste0("flowcytometer_traits_bacteria.", plate, ".rds")))
    )
  }
  saveRDS(
    traits,
    file = file.path(add_path, paste0("flowcytometer_traits_bacteria.rds"))
  )
  utils::write.csv(
    traits,
    file = file.path(add_path, paste0("flowcytometer_traits_algae.csv")),
    row.names = FALSE
  )
  #
  traits <- NULL
  #
  for (plate in plates){
    traits <- rbind(
      traits,
      readRDS(file.path(add_path, paste0("flowcytometer_traits_algae.", plate, ".rds")))
    )
  }
  saveRDS(
    traits,
    file = file.path(add_path, paste0("flowcytometer_traits_algae.rds"))
  )
  utils::write.csv(
    traits,
    file = file.path(add_path, paste0("flowcytometer_traits_algae.csv")),
    row.names = FALSE
  )


  # Finalize ----------------------------------------------------------------

  unlink(processing)
  message("   done")
  message("########################################################")

  invisible(TRUE)
}
