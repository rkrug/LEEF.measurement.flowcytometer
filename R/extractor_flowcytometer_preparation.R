#' Extractor flowcytometer preparation data
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
  message("Extracting flowcytometer preparation ...\n")

  ##
  suppressWarnings(
    {
      processing <- file.path(normalizePath(output), "flowcytometer", paste0("EXTRACTING.FLOWCYTOMETER.PREPARATION", ".PROCESSING"))
      error <- file.path(normalizePath(output), "flowcytometer", paste0("ERROR.EXTRACTING.FLOWCYTOMETER.PREPARATION", ".ERROR"))
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

  # Based on flowcyt_1_c6_to_RData.R ----------------------------------------
  # Converting the Flowcytometer Output of bacterial abundances into a usable data frame
  # David Inauen, 19.06.2017


  # Get fcs file names ------------------------------------------------------

  fcs_path <- file.path( input, "flowcytometer" )
  fcs_path <- list.dirs(
    fcs_path,
    full.names = TRUE,
    recursive = FALSE
  )
  fcs_files <- list.files(
    path = fcs_path,
    pattern = "*.fcs",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(fcs_path) != 1) {
    unlink(processing)
    message("Exactly one directory is expected in the flowcytometer folder\n")
    message("\n########################################################\n")
    return(invisible(FALSE))
  }

  if (length(fcs_files) == 0) {
    unlink(processing)
    message("nothing to extract\n")
    message("\n########################################################\n")
    return(invisible(TRUE))
  }


  # Load parameter files\ ----------------------------------------------------

  metadata <- read.csv(file.path( input, "flowcytometer", "metadata_flowcytometer.csv" ))

  # check file sizes and delete empty wells ---------------------------------

  fcs_files <- fcs_files[ file.info(fcs_files)[,"size"] > 3000 ]

  # read flowSet automatically ----------------------------------------------

  # Due to changes in flowCore (simplifications) the code should be changed.
  # See https://support.bioconductor.org/p/p132747/#p132763 for details

  # Old Code
  # fsa <- flowCore::read.flowSet(
  #   path = fcs_path,
  #   transformation = FALSE,
  #   phenoData = list(
  #     filename = "#SAMPLE",
  #     sample = "$SMNO",
  #     date = "$DATE",
  #     volume = "$VOL",
  #     proj = "$PROJ"
  #   )
  # )

  ### begin New code from link above
  # read fcs
  fsa <- flowCore::read.flowSet(
    path = fcs_path
  )

  #extract keyword list
  kw <- flowCore::keyword(
    fsa,
    keyword = list(
      # filename = "#SAMPLE",
      sample = "$WELLID", ## "$SMNO",
      # date = "$DATE",
      volume = "$VOL",
      proj = "$PROJ" ## needed?????
    )
  )

  #assign it to pdata of fs
  flowCore::pData(fsa) <- as.data.frame(kw)


  # CREATE FLOW DATA FRAME AND FILL WITH UNGATED COUNT ----------------------

  flow.data <- flowCore::pData( flowCore::phenoData(fsa) )

  # find the number of events (equals the number of rows)
  num <- sapply(
    1:length(fsa),
    function(i) {
      num <- dim( flowCore::exprs(fsa[[i]]) )[1]
    }
  )
  flow.data <- cbind(
    flow.data,
    "total.counts" = num
  )

  # Extract the volume of medium sampled
  flow.data[["volume"]] <- as.numeric(
    as.character(
      flowCore::phenoData(fsa)[["volume"]]
    )
  )

  # RL: new: merge with metadata
  flow.data <- merge(flow.data, metadata, by = "sample", all.x = TRUE)


  # calculate events recorded per ml
  #RL: new: use dilution_factor from metadata
  flow.data[["tot_density_perml"]] <- flow.data[["total.counts"]] * 1000000 / flow.data[["volume"]] * flow.data[["dilution_factor"]]
  #   flow.data[["specname"]] <- paste(
  #     flow.data[["filename"]],
  #     flow.data[["proj"]],
  #     sep = "_"
  #   )


  # standardize naming
  # flow.data <- flow.data[, c("filename","sample","date","volume","total.counts","tot_density_perml","specname")]
  flow.data <- flow.data[, c("sample", "bottle", "volume","total.counts","tot_density_perml", "dilution_factor")]

  rownames(flow.data) <- NULL


  # define transformation
  # logTrans <- logTransform( transformationId="log10-transformation", logbase = 10, r = 1, d = 1 )
  # aTrans <- truncateTransform( "truncate at 1", a = 1 )

  # exclude values < 1 (RL: use FL1-A and FL3-A instead of -H; added line for FL4-A)
  fsa <- flowCore::transform(
    fsa,
    `FL1-A` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`FL1-A`),
    `FL3-A` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`FL3-A`),
    `FL4-A` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`FL4-A`),
    `FSC-A` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`FSC-A`),
    `SSC-A` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`SSC-A`),
    `Width` = flowCore::truncateTransform( "truncate at 1", a = 1 )(`Width`)
  )
  # log transform (RL: use FL1-A and FL3-A instead of -H; added line for FL4-A)
  fsa <- flowCore::transform(
    fsa,
    `FL1-A` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`FL1-A`),
    `FL3-A` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`FL3-A`),
    `FL4-A` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`FL4-A`),
    `FSC-A` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`FSC-A`),
    `SSC-A` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`SSC-A`),
    `Width` = flowCore::logTransform( transformationId = "log10-transformation", logbase = 10, r = 1, d = 1 )(`Width`)
  )


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
    file = file.path(add_path, "flowcytometer_ungated.csv"),
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
