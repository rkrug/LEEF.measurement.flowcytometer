#' Extractor flowcytometer data
#'
#'
#' This function is extracting data to be added to the database
#' (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param particles \code{character} vector containing the groups to extract.
#'   Supported are \code{"bacteria"}, \code{"LNA"}, \code{"MNA"}, \code{"HNA"},
#'   \code{"algae"} and \code{"all"} which does no gating.
#' @param metadata_flowcytometer the content of the file \code{metadata_flowcytometer.csv} which will be linked into the traits
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
#' @export
#'
extract_traits <- function(
    input,
    particles = c("bacteria", "LNA", "MNA", "HNA", "algae"),
    metadata_flowcytometer
) {

  # function to gate each plate ---------------------------------------------

  traits <- function(
    input,
    particles
  )  {
    # extraction function -----------------------------------------------------

    extr_traits <- function(
    pop
    ){
      traits <- flowCore::fsApply(
        pop,
        function(p){
          result <- list()
          result$sample <- unlist(flowCore::keyword(p, "$WELLID"))
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
      traits$timestamp <- yaml::read_yaml(file.path(input, "sample_metadata.yml"))$timestamp
      return(traits)
    }

    #  read data files ----------------------------------------------------------

    message("   reading data  ...")
    fsa <- readRDS(file = file.path(input, paste0("flowcytometer_fsa_ungated.rds")))

    gates <- calculate_gates(input)

    # The Subsetting and Extraction ------------------------------------------------------------

    result <- list()

    bacteria_pop <- NULL

    if ("all" %in% particles){
      result$all <- extr_traits(fsa)
    }
    if ("bacteria" %in% particles){
      message("   gating bacteria ...")
      bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)
      result$bacteria <- extr_traits(bacteria_pop)
    }
    if ("LNA" %in% particles){
      message("   gating LNA ...")
      if (is.null(bacteria_pop)){
        bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)
      }
      result$LNA <- extr_traits(Subset(bacteria_pop, gates$bacteria$rg_LNA))
    }
    if ("MNA" %in% particles){
      message("   gating MNA ...")
      if (is.null(bacteria_pop)){
        bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)
      }
      result$MNA <- extr_traits(Subset(bacteria_pop, gates$bacteria$rg_MNA))
    }
    if ("HNA" %in% particles){
      message("   gating HNA ...")
      if (is.null(bacteria_pop)){
        bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)
      }
      result$HNA <- extr_traits(Subset(bacteria_pop, gates$bacteria$rg_HNA))
    }
    if ("algae" %in% particles){
      message("   gating algae ...")
      result$algae <- extr_traits(Subset(fsa,  gates$algae$algae_gate))
    }

    result <- lapply(
      result,
      function(traits){
        traits <- merge(
          traits,
          metadata_flowcytometer,
          by.x = c("sample"),
          by.y = c("sample"),
          all.x = TRUE,
          all.y = FALSE
        )
        return(traits)
      }
    )

    return(result)
  }


  # Do the trait extraction --------------------------------------------

  results <- traits(input, particles)

  # Finalize ----------------------------------------------------------------

  message("   done")
  message("########################################################")

  return(results)
}
