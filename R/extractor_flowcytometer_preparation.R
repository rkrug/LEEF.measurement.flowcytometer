#' Extractor flowcytometer preparation data
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
#' @importFrom flowCore read.flowSet pData phenoData exprs logTransform truncateTransform transform rectangleGate Subset transformList
#' @importFrom  yaml read_yaml
#' @importFrom  stats setNames
#' @importFrom utils read.csv write.csv
#' @export
#'
extractor_flowcytometer_preparation <- function(
    input,
    output
) {
  add_path <- file.path(output, "flowcytometer")
  dir.create(add_path, recursive = TRUE, showWarnings = FALSE)
  log.file <- tempfile(pattern="flowcytometer.", tmpdir = add_path, fileext = ".log")
  loggit::set_logfile(log.file)

  message("########################################################")
  message("   preparing flowcytometer...")

  ##
  suppressWarnings({
    processing <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("EXTRACTING.FLOWCYTOMETER.PREPARATION", ".PROCESSING")
    )
    error <- file.path(
      normalizePath(output), "flowcytometer",
      paste0("ERROR.EXTRACTING.FLOWCYTOMETER.PREPARATION", ".ERROR")
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


  # Function to extract individual fcs folder per plate / folder -------------------------------


  extract_plate <- function(
    plate,
    input,
    output
  ){
    # Based on flowcyt_1_c6_to_RData.R ----------------------------------------
    # Converting the Flowcytometer Output of bacterial abundances into a usable data frame
    # David Inauen, 19.06.2017
    #
    # Get fcs file names ------------------------------------------------------
    fcs_path <- file.path(input, "flowcytometer", plate)

    fcs_files <- list.files(
      path = fcs_path,
      pattern = "*.fcs",
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(fcs_files) == 0) {
      unlink(processing)
      message("   nothing to extract")
      message("########################################################")
      return(invisible(TRUE))
    }


    # Load parameter files ----------------------------------------------------

    metadata <- utils::read.csv(file.path(input, "flowcytometer", "metadata_flowcytometer.csv"))
    metadata <- metadata[metadata$plate == plate,]

    #############################################################
    # <<<< BEGIN SCRIPT   #######################################
    #############################################################

    # check file sizes and exclude empty wells ---------------------------------

    fcs_files <- fcs_files[file.size(fcs_files) > 4300]

    # read flowSet automatically ----------------------------------------------

    # read fcs
    fsa <- flowCore::read.flowSet(
      files = fcs_files
    )


    #  BEGIN from script --------------------------------------------------------------------

    # extract keyword list
    # Uriah: I put date back in because it is needed in the script below...
    # or can I remove it in the script below?
    kw <- flowCore::keyword(
      fsa,
      keyword = list(
        sample = "$WELLID",
        date = "$DATE",
        volume = "$VOL",
        proj = "$PROJ"
      )
    ) ## PROJ needed?????

    ###### ADDED FOR EAWAG MEASUREMENTS BEGIN ######

    if (!("sample" %in% colnames(kw))){
      sample <- rownames(kw)
      sample <- gsub(".fcs$", "", sample)
      kw <- cbind(sample = sample, kw)
    }

    ###### ADDED FOR EAWAG MEASUREMENTS END ######

    kw <- cbind(filename = rownames(kw), kw)

    #assign it to pdata of fs
    flowCore::pData(fsa) <- as.data.frame(kw)

    # CREATE FLOW DATA FRAME AND FILL WITH UNGATED COUNT ----------------------
    flow.data <- flowCore::pData(flowCore::phenoData(fsa))

    # find the number of events (equals the number of rows)
    num <- sapply(
      seq_along(fsa),
      function(i) {
        dim(flowCore::exprs(fsa[[i]]))[1]
      }
    )
    flow.data <- cbind(flow.data, "total.counts" = num)

    # Extract the volume sampled
    flow.data$volume <- as.numeric(as.character(flowCore::phenoData(fsa)$volume))

    # RL: new: merge with metadata
    flow.data <- merge(flow.data, metadata, )

    # calculate events recorded per ml
    #RL: new: use dilution_factor from metadata
    flow.data$tot_density_perml <- flow.data$total.counts * 1000000 / flow.data$volume * flow.data$dilution_factor
    flow.data$specname <- paste(flow.data$filename, flow.data$proj, sep = "_")

    flow.data <- flow.data[, c("filename", "bottle", "date", "sample", "volume", "total.counts",
                               "tot_density_perml", "specname", "dilution_factor")]
    rownames(flow.data) <- NULL

    # exclude values < 1 (RL: use FL1-A and FL3-A instead of -H; added line for FL4-A)
    fsa <- flowCore::transform(
      fsa,
      flowCore::transformList(
        c("FL1-A", "FL3-A", "FL4-A", "FSC-A", "SSC-A", "Width"),
        truncateTransform("truncate at 1")
      )
    )

    # log transform (RL: use FL1-A and FL3-A instead of -H; added line for FL4-A)
    fsa <- flowCore::transform(
      fsa,
      flowCore::transformList(c("FL1-A", "FL3-A", "FL4-A", "FSC-A", "SSC-A", "Width"), "log10")
    )

    #############################################################
    # >>>> END SCRIPT   #########################################
    #############################################################

    # SAVE --------------------------------------------------------------------

    add_path <- file.path(output, "flowcytometer")
    dir.create(add_path, recursive = TRUE, showWarnings = FALSE)

    timestamp <- yaml::read_yaml(file.path(input,  "flowcytometer", "sample_metadata.yml"))$timestamp

    flow.data <- cbind(
      timestamp = timestamp,
      plate = plate,
      flow.data
    )

    saveRDS(
      flow.data,
      file = file.path(output, "flowcytometer", paste0("flowcytometer_ungated.", plate, ".rds"))
    )

    saveRDS(
      fsa,
      file = file.path(output, "flowcytometer", paste0("flowcytometer_fsa_ungated.", plate, ".rds"))
    )
  }


# Identify all fcs folder and iterate through them ------------------------


  plates <- grep(
    "p_",
    list.dirs(
      file.path(input, "flowcytometer"),
      full.names = FALSE
    ),
    value = TRUE
  )

  lapply(
    plates,
    extract_plate,
    input = input,
    output = output
  )


# Do final copying of `.ciplus` files -------------------------------------


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
