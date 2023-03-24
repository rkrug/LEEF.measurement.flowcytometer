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
#' @param metadata_flowcytometer the content of the file
#'   \code{metadata_flowcytometer.csv} which will be linked into the traits
#' @param excl_FSCA_0 boolean. If \code{TRUE}, \code{FSA.A <= 0} will be fitered
#'   out by using a rectangular filter
#'   \code{flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" =
#'   c(0.00000000001, +Inf))}
#' @param use_H if \code{TRUE}, gating will be done using \code{height},
#'   otherwie \code{area}
#' @param timestamp timestamp. Default: read from \code{sample_metadata.yml}
#' @param fsa if \code{NULL}, \code{fsa} will be read in, otherwise the
#'   \code{fsa}
#' @param gates_coordinates if \code{NULL}, \code{gates_coordinates} will be
#'   read in, otherwise the \code{gates_coordinates}
#' @param wellid_keyword the kwyword which is used to identify the well ID.
#'   Usually "$WELLID" (default), but for the EAWAG Flowcytometer it is "$SMNO".
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
    input = NULL,
    particles = c("bacteria", "LNA", "MNA", "HNA", "algae"),
    metadata_flowcytometer,
    excl_FSCA_0 = FALSE,
    use_H = FALSE,
    timestamp = yaml::read_yaml(file.path(input, "sample_metadata.yml"))$timestamp,
    fsa = NULL,
    gates_coordinates = NULL,
    wellid_keyword = "$WELLID"
) {

  # function to gate each plate ---------------------------------------------

  traits <- function(
    input,
    particles,
    wellid_keyword
  )  {
    # extraction function -----------------------------------------------------

    extr_traits <- function(
      pop,
      wellid_keyword
    ){
      traits <- flowCore::fsApply(
        pop,
        function(p){
          result <- list()
          result$sample <- unlist(flowCore::keyword(p, wellid_keyword))
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
      traits$timestamp <- timestamp
      return(traits)
    }

    #  read data files ----------------------------------------------------------

    message("   reading data  ...")
    if (is.null(fsa)){
      fsa <- readRDS(file = file.path(input, paste0("flowcytometer_fsa_ungated.rds")))
    }

    if (excl_FSCA_0){
      g0 <- flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(0.00000000001, +Inf))
      fsa <- flowCore::Subset(fsa, g0)
    }

    if (is.null(gates_coordinates)){
      gates_coordinates <- utils::read.csv(file.path(input_dir, "gates_coordinates.csv"))
    }

    if (use_H) {
      gates <- calculate_gates_H(gates_coordinates = gates_coordinates)
    } else {
      gates <- calculate_gates(gates_coordinates = gates_coordinates)
    }

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
      result$algae <- extr_traits(Subset(fsa,  gates$algae$algae$algae_gate))
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

  results <- traits(input, particles, wellid_keyword)

  # Finalize ----------------------------------------------------------------

  message("   done")
  message("########################################################")

  return(results)
}
