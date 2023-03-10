#' Calculate gates based on gate co-ordinates in csv file using \code{area}
#'
#' @param input_dir directory which contains the file \code{gates_coordinates.csv}
#' @param gates_coordinates contains gate co-ordinates. If \code{NULL} *default)
#'   the gates cootdinates are read in from \code{input_dir}
#'
#' @return \code{list} containing two objects: gates for bacteria and gates for algae
#'
#' @importFrom utils read.csv
#' @importFrom  flowCore polygonGate rectangleGate
#' @export
#'
#' @examples
calculate_gates_H <- function(
    input_dir = NULL,
    gates_coordinates = NULL
){
  if (is.null(gates_coordinates)) {
    gates_coordinates <- utils::read.csv(file.path(input_dir, "gates_coordinates.csv"))
  }

  gates <- list()

  # bacteria gate
  polyGate_bacteria <- as.matrix(gates_coordinates[1:4, 1:2])
  colnames(polyGate_bacteria) <- c("FL1-H", "FL3-H")
  bacteria_gate <- flowCore::polygonGate(filterId = "Bacteria", .gate = polyGate_bacteria)

  # gate for different size classes of bacteria
  LNA_coordinates <- as.matrix(gates_coordinates[, 3])
  LNA_coordinates <- na.omit(LNA_coordinates)
  colnames(LNA_coordinates) <- c("FL1-H")

  MNA_coordinates <- as.matrix(gates_coordinates[, 4])
  MNA_coordinates <- na.omit(MNA_coordinates)
  colnames(MNA_coordinates) <- c("FL1-H")

  HNA_coordinates <- as.matrix(gates_coordinates[, 5])
  HNA_coordinates <- na.omit(HNA_coordinates)
  colnames(HNA_coordinates) <- c("FL1-H")

  rg_LNA <- flowCore::rectangleGate("FL1-H" = LNA_coordinates, filterId = "LNA")
  rg_MNA <- flowCore::rectangleGate("FL1-H" = MNA_coordinates, filterId = "MNA")
  rg_HNA <- flowCore::rectangleGate("FL1-H" = HNA_coordinates, filterId = "HNA")

  gates$bacteria <- list(
    bacteria_gate = bacteria_gate,
    rg_LNA = rg_LNA,
    rg_MNA = rg_MNA,
    rg_HNA = rg_HNA
  )

  # algae gate
  polyGate_algae <- as.matrix(gates_coordinates[1:4, 6:7])
  colnames(polyGate_algae) <- c("FL1-H", "FL4-H")
  algae_gate <- flowCore::polygonGate(filterId = "Algae", .gate = polyGate_algae)

  gates$algae <- list(
    algae_gate = algae_gate
  )

  return(gates)
}

