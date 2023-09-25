#' Extract flowcytometer traits
#'
#' NOT USED IN PIPLELINE!!!!
#' This function is extracting data to be added to the database
#' (and therefore make accessible for further analysis and forecasting)
#' from \code{.fcs} files.
#'
#' @param input directory from which to read the data
#' @param metadata_flowcytometer the content of the file
#'   \code{metadata_flowcytometer.csv} which will be linked into the traits
#' @param min_FSC.A numeric. If \code{!NULL}, \code{FSA.A <= min_FSC.A} will be fitered out by using
#'   a rectangular filter
#'   \code{flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))}
#' @param use_H if \code{TRUE}, gating will be done using \code{height},
#'   otherwie \code{area}
#' @param timestamp timestamp. Default: read from \code{sample_metadata.yml}
#' @param fsa_p1 if \code{NULL}, \code{fsa_p1} will be read in, otherwise the
#'   \code{fsa}
#' @param fsa_p2 if \code{NULL}, \code{fsa_p2} will be read in, otherwise the
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
#' @import loggit
#' @export
#'
extract_traits <- function(
    input = NULL,
    particles = c("bacteria", "LNA", "MNA", "HNA", "algae"),
    metadata_flowcytometer,
    min_FSC.A = NULL,
    use_H = FALSE,
    timestamp = yaml::read_yaml(file.path(input, "sample_metadata.yml"))$timestamp,
    fsa_p1 = NULL,
    fsa_p2 = NULL,
    gates_coordinates = NULL,
    wellid_keyword = "$WELLID") {
    message("   reading data  ...")

    if (is.null(fsa_p1)) {
        fsa_p1 <- readRDS(file = file.path(input, paste0("flowcytometer_fsa_ungated.p_1.rds")))
    }

    if (is.null(fsa_p2)) {
        fsa_p2 <- readRDS(file = file.path(input, paste0("flowcytometer_fsa_ungated.p_1.rds")))
    }

    if (!is.null(min_FSC.A)) {
        g0 <- flowCore::rectangleGate(filterId = "filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))
        fsa_p1 <- flowCore::Subset(fsa_p1, g0)
        fsa_p2 <- flowCore::Subset(fsa_p2, g0)
    }

    if (is.null(gates_coordinates)) {
        gates_coordinates <- utils::read.csv(file.path(input, "gates_coordinates.csv"))
    }

    if (use_H) {
        gates <- calculate_gates_H(gates_coordinates = gates_coordinates)
    } else {
        gates <- calculate_gates(gates_coordinates = gates_coordinates)
    }



    # function to gate each plate ---------------------------------------------

    trait_plate <- function(plate_id,
                            fsa,
                            gates = gates,
                            wellid_keyword = wellid_keyword,
                            timestamp = timestamp) {
        # extraction function -----------------------------------------------------

        extr_traits <- function(pop,
                                wellid_keyword = wellid_keyword,
                                timestamp = timestamp) {
            traits <- flowCore::fsApply(
                pop,
                function(p) {
                    result <- list()
                    result$sample <- unlist(flowCore::keyword(p, wellid_keyword))
                    result$plate <- plate_id
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

        traits <- list()
        # The Bacteria ------------------------------------------------------------

        bacteria_pop <- Subset(fsa, gates$bacteria$bacteria_gate)

        # LNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_LNA)
        # MNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_MNA)
        # HNA_pop <- Subset(bacteria_pop, gates$bacteria$rg_HNA)

        traits$bacteria <- extr_traits(bacteria_pop, wellid_keyword = wellid_keyword, timestamp = timestamp)

        # The Algae ---------------------------------------------------------------

        algae_pop <- Subset(fsa, gates$algae$algae_gate)

        traits$algae <- extr_traits(bacteria_pop, wellid_keyword = wellid_keyword, timestamp = timestamp)

        return(traits)
    }


    # Do the trait extraction for all plates --------------------------------------------

    p1 <- trait_plate(
        plate_id = "p_1",
        fsa = fsa_p1,
        gates = gates,
        timestamp = timestamp,
        wellid_keyword = wellid_keyword
    )
    p2 <- trait_plate(
        plate_id = "p_2",
        fsa = fsa_p2,
        gates = gates,
        timestamp = timestamp,
        wellid_keyword = wellid_keyword
    )

    traits <- list()

    traits$bacteria <- rbind(
        p1$bacteria,
        p2$bacteria
    )
    traits$algae <- rbind(
        p1$algae,
        p2$algae
    )

    traits$bacteria <- merge(
        traits$bacteria,
        metadata_flowcytometer,
        by.x = c("sample", "plate"),
        by.y = c("sample", "plate"),
        all.x = TRUE,
        all.y = FALSE
    )
    traits$algae <- merge(
        traits$algae,
        metadata_flowcytometer,
        by.x = c("sample", "plate"),
        by.y = c("sample", "plate"),
        all.x = TRUE,
        all.y = FALSE
    )



    # Finalize ----------------------------------------------------------------

    message("   done")
    message("########################################################")

    return(traits)
}
