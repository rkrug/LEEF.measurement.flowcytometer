#' Calculate gates based on gate co-ordinates in csv file
#'
#' @param input_dir directory which contains the file \code{gates_coordinates.csv}
#'
#' @return \code{list} containing two objects: gates for bacteria and gates for algae
#'
#' @importFrom utils read.csv
#' @importFrom  flowCore polygonGate rectangleGate
#' @export
#'
#' @examples
calculate_gates <- function(
    input_dir
){
  gates_coordinates <- utils::read.csv(file.path(input_dir, "gates_coordinates.csv"))

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

  return(gates)
}

