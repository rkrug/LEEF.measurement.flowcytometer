#' Get density using \code{area} for gating
#'
#' @param gates_coordinates gates co-ordinates as read from the file \code{gates_coordinates.csv}
#' @param fsa as read from \code{flowcytometer_fsa_ungated.rds}
#' @param flow.data as read from \code{flowcytometer_ungated.csv}
#' @param min_FSC.A numeric. If \code{!NULL}, \code{FSA.A <= min_FSC.A} will be fitered out by using
#'   a rectangular filter
#'   \code{flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))}
#'
#' @return  gated \code{flow.data} density as saved in \code{flowcytometer_density.csv} and actual gates used
#' @export
#'
#' @examples
dens <- function(
    gates_coordinates,
    fsa,
    flow.data,
    min_FSC.A = 0.0000000001
){

  #RL: defining gates

  gates <- calculate_gates( gates_coordinates = gates_coordinates)


  # #----- HAVING A LOOK AT THE GATING -----#

  # i <- 1
  # flowViz::xyplot(`FL3-A` ~ `FL1-A`, data = fsa[[i]], filter = gates$bacteria$bacteria_gate)
  # # fsa[[i]]
  # flowViz::densityplot(~ `FL1-A`, data = fsa[[i]], filter = rg_LNA)
  # # fsa[[i]]
  # flowViz::xyplot(`FL4-A` ~ `FL1-A`, data = fsa[[i]], filter = algae_gate)
  # # fsa[[i]]

  if (!is.null(min_FSC.A)){
    g0 <- flowCore::rectangleGate(filterId="filter_out_0", "FSC-A" = c(min_FSC.A, +Inf))
    fsa <- flowCore::Subset(fsa, g0)
  }
  
  # ABUNDANCE DYNAMICS ------------------------------------------------------

  # applying filter to whole flowSet

  result <- flowCore::filter(fsa, gates$bacteria$bacteria_gate)

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
  subset.bacteria <- flowCore::Subset(fsa, gates$bacteria$bacteria_gate)

  # applying filter to bacteria to get the three bacteria populations
  LNA <- flowCore::filter(subset.bacteria, gates$bacteria$rg_LNA)
  MNA <- flowCore::filter(subset.bacteria, gates$bacteria$rg_MNA)
  HNA <- flowCore::filter(subset.bacteria, gates$bacteria$rg_HNA)

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
  algae <- flowCore::filter(
    flowCore::Subset(fsa, !gates$bacteria$bacteria_gate),
    gates$algae$algae_gate
  )

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


  return(
    list(
      flow.data = flow.data,
      gates = gates
      )
    )
}
