#define globals to suppress github actions warnings
utils::globalVariables(
  names = c('x_global_px', 'y_global_px', 'cellID', 'cell_ID', 'fov'),
  package = 'Flat2OMETiff',
  add = TRUE
)

#' @title Flat2Matrix
#'
#' @description This function converts flat files from CosMx into a directory of sparse matrices. One file for each channel will be created within the `outputDirectory` directory.
#'
#' @param expressionCsvPath Path to the expression matrix CSV file.
#' @param polygonCsvPath Path to the polygon CSV file.
#' @param outputDirectory Path to the directory where the sparse matrices will be written.
#' @param channelNamesAreExpressionCSVColumnNames Boolean to control whether the channel names are taken from the expression CSV file column names. If not, the channels are indexed numerically.
#' @importFrom data.table fread
#' @importFrom dplyr filter select mutate
#'
#' @export

#TODO: parallelize
Flat2Matrix <- function(expressionCsvPath = NULL,
                        polygonCsvPath = NULL,
                        outputDirectory = "./sparse_matrix_files/",
                        channelNamesAreExpressionCSVColumnNames = TRUE) {
  #if the base path doesn't exist, create it
  if (!dir.exists(outputDirectory)) {
    dir.create(outputDirectory)
  }

  #read in the flat files
  flat_file_list <- .ReadFlatFiles(expressionCsvPath = expressionCsvPath,
                                  polygonCsvPath = polygonCsvPath)
  flat_file_list$polygons <- flat_file_list$polygons %>%
    mutate(x_global_px = round(x_global_px),
           y_global_px = round(y_global_px))


  #populate image bounds
  x_min <- min(flat_file_list$polygons$x_global_px)
  x_max <- max(flat_file_list$polygons$x_global_px)
  y_min <- min(flat_file_list$polygons$y_global_px)
  y_max <- max(flat_file_list$polygons$y_global_px)

  #offset the image if it has negative coordinates for some CosMx-internal reason.
  if (x_min < 0 | y_min < 0) {
    offset_y <- ifelse(y_min < 0, abs(y_min) + 1, 0)
    offset_x <- ifelse(x_min < 0, abs(x_min) + 1, 0)
    x_min <- x_min + offset_x
    x_max <- x_max + offset_x
    y_min <- y_min + offset_y
    y_max <- y_max + offset_y
  } else {
    offset_x <- 0
    offset_y <- 0
  }

  #harvest protein names from the expression matrix, but drop the cellular/FOV identifiers
  channel_names <- unique(colnames(flat_file_list$expression))[!(unique(colnames(flat_file_list$expression)) %in% c("fov", "cell_ID"))]

  #TODO: remove the 1:3 when the function is ready for production
  for (channel in channel_names[1:3]) {
    #write blank image with correct global dimensions
    sparseMatrix <- .InitalizeSparseMatrix(x_min, x_max, y_min, y_max)

    #this sets up a vector of cells within FOVs that can be strsplit via:
    #FOV ID = strsplit(unique_cell_ids, split ="_")[2]
    #Cell ID = strsplit(unique_cell_ids, split ="_")[4]
    unique_cells <- paste0("fov_", flat_file_list$expression$fov, "_cellID_", flat_file_list$expression$cell_ID)

    #iterate over cells and fill the convex hulls (defined by the polygon flat files) with expression data
    for (cell in unique_cells){
      fov_ID <-  unlist(strsplit(cell, split ="_"))[2]
      cell_index <- unlist(strsplit(cell, split ="_"))[4]

      vertices <- flat_file_list$polygons %>%
          dplyr::filter(cellID == cell_index & fov == fov_ID) %>%
          dplyr::select(x_global_px, y_global_px) %>%
          dplyr::mutate(x_global_px = x_global_px + offset_x,
                        y_global_px = y_global_px + offset_y) %>%
          as.matrix()

      #create hull
      hull <- geometry::convhulln(cbind(vertices[,1], vertices[,2]), options = "FA")

      #get a bounding box for the cell
      cell_x_min <- min(vertices[,1])
      cell_x_max <- max(vertices[,1])
      cell_y_min <- min(vertices[,2])
      cell_y_max <- max(vertices[,2])

      #create a grid of points within the cell's bounding box
      cell_sub_image <- expand.grid(x = seq(cell_x_min, cell_x_max, 1), y = seq(cell_y_min, cell_y_max, 1))

      #look up the intensity of the channel in the cell's expression matrix
      #TODO / NOTE: This will eventually probably scale to the transcriptomics flat files,
      #in which this function won't have a single value per cell, but that data will be very dense,
      #so I'll need to know that all of this works and is use before taking that on.
      cell_expression <- flat_file_list$expression %>%
        dplyr::filter(fov == fov_ID & cell_ID == cell_index) %>%
        dplyr::select(dplyr::all_of(channel)) %>%
        as.numeric()
      #define the pixels inside of and on the hull
      inside <- sp::point.in.polygon(cell_sub_image[,1], cell_sub_image[,2], hull$p[,1], hull$p[,2])
      #write the expression of the current protein/gene as the intensity of the interior of the cell
      #TODO: this will need to be adjusted if we move from average intensities to per-pixel/event data.
      inside[inside == 1] <- cell_expression

      #if it's a border/vertex pixel, set the maximum intensity so that the border is visible across every channel
      inside <- ifelse( inside == 2 | inside == 3, yes = 255, no = inside)

      #define the sub-image in units of pixels in the larger image
      i_indices <- cell_sub_image$x
      j_indices <- cell_sub_image$y

      #filter out the zero values to keep the matrix sparse, but in coordinates of the sub-image
      non_zero_indices <- which(inside != 0)

      #filter out indices in the dense sub-image that are zero.
      global_i_indices <- i_indices[non_zero_indices]
      global_j_indices <- j_indices[non_zero_indices]

      #check for and handle zero values
      #the zeroes in the sparse matrix will yield NA, since they're still sparse, so . != 0 yields NA
      valid_indices <- !is.na(global_i_indices) & !is.na(global_j_indices)


      #filter valid indices
      valid_global_i_indices <- global_i_indices[valid_indices]
      valid_global_j_indices <- global_j_indices[valid_indices]
      valid_inside_values <- inside[non_zero_indices][valid_indices]
      #the indices of the sub image aren't guaranteed to start at 0/1, but writing to the sparse matrix requires that the indices fit into:
      # 1:x_max-x_min and 1:y_max-y_min, so we need to adjust the indices to fit into those ranges.
      valid_i_adjusted <- valid_global_i_indices - x_min
      valid_j_adjusted <- valid_global_j_indices - y_min


      #convert the cell's dense sub-image into a sparse matrix with units of the original image for matrix addition.
      cell_sparse_matrix <- Matrix::sparseMatrix(
        i = valid_i_adjusted,
        j = valid_j_adjusted,
        x = valid_inside_values,
        dims = c(length(seq(x_min, x_max, 1)),
                 length(seq(y_min, y_max, 1))),
        dimnames = list(seq(x_min, x_max, 1) ,
                 seq(y_min, y_max, 1))
      )
      #matrix addition to write the cell's intensity information into the initialized image for the channel
      sparseMatrix[valid_i_adjusted, valid_j_adjusted] <-
        sparseMatrix[valid_i_adjusted, valid_j_adjusted] +
        cell_sparse_matrix[valid_i_adjusted, valid_j_adjusted]
    }
    if (channelNamesAreExpressionCSVColumnNames) {
      channel_name <- channel
    } else {
      channel_name <- which(channel_names == channel)
    }
    #write the channel's sparse matrix to disk
    Matrix::writeMM(sparseMatrix, file = paste0(outputDirectory, "/sparse_matrix_", channel_name, ".mtx"))
  }
}

#' @title InitalizeSparseMatrix
#' @description This function initializes a sparse matrix with the dimensions of the global image.
#' @param x_min The minimum x-coordinate of the global image.
#' @param x_max The maximum x-coordinate of the global image.
#' @param y_min The minimum y-coordinate of the global image.
#' @param y_max The maximum y-coordinate of the global image.
#' @return A sparse matrix with the dimensions of the CosMx image.

.InitalizeSparseMatrix <- function(x_min, x_max, y_min, y_max) {
  sparse_matrix <- Matrix::Matrix(0,
                 nrow = length(seq(x_min, x_max, 1)),
                 ncol = length(seq(y_min, y_max, 1)),
                 dimnames = list(seq(x_min, x_max, 1), seq(y_min, y_max, 1)),
                 sparse = TRUE)
  return(sparse_matrix)
}

#' @title ReadFlatFiles
#' @description This function reads in the flat files from CosMx.
#' @param expressionCsvPath Path to the expression matrix CSV file.
#' @param polygonCsvPath Path to the polygon CSV file.
#' @return A list containing the expression and polygon data tables.

.ReadFlatFiles <- function(expressionCsvPath = NULL,
                          polygonCsvPath = NULL) {
  expression <- data.table::fread(expressionCsvPath)
  polygons <- data.table::fread(polygonCsvPath)
  return(list(expression = expression, polygons = polygons))
}

#' @title WriteTiffChannels
#' @description This function writes the sparse matrices to tiff files using the python tifffile library.
#' @param sparseMatrixDirectory Path to the directory containing the sparse matrices.
#' @param outputDirectory Path to the directory where the tiff files will be written.
#' @export

WriteTiffChannels <- function(sparseMatrixDirectory = "./sparse_matrix_files/",
                              outputDirectory = "./tiff_channels/") {
  #enforce absolute paths
  sparseMatrixDirectory <- R.utils::getAbsolutePath(sparseMatrixDirectory)
  outputDirectory <- R.utils::getAbsolutePath(outputDirectory)
  for (matrixFile in list.files(sparseMatrixDirectory, full.names = T)) {
    print(matrixFile)
    start_time <- Sys.time()
    script_contents <- readr::read_file(system.file("scripts/WriteChannels.py", package = "Flat2OMETiff"))
    script <- tempfile()

    #write the function definition to the temp file
    readr::write_file(script_contents, script)
    function_call <- paste0("WriteChannelTiff(sparseMatrixFile = '", matrixFile,
                     "',outputDirectory = '", outputDirectory,
                     "')")
    print(function_call)
    #write the function call with arguments to the end of the script and execute
    readr::write_file(function_call, script, append = TRUE)
    system2(reticulate::py_exe(), script)
    end_time <- Sys.time()
    print(paste0("Time to write ", matrixFile, " to tiff: ", end_time - start_time))
  }


}

#' @title StitchTiffChannels
#' @description This function stitches the tiff channels together using the python tifffile library.
#' @param tiffChannelDirectory Path to the directory containing the tiff channels.
#' @param outputDirectory Path to the directory where the stitched tiff files will be written.
#' @param channelWhitelist A vector of channel names to include in the stitched tiff files.
#' @export

StitchTiffChannels <- function(tiffChannelDirectory = "./tiff_channels/",
                               outputDirectory = "./stitched_tiff/",
                              channelWhitelist = NULL) {
  #TODO: add a check for the channelWhitelist to ensure that it's a vector and create a default behavior.
  #enforce absolute paths
  tiffChannelDirectory <- R.utils::getAbsolutePath(tiffChannelDirectory)
  outputDirectory <- R.utils::getAbsolutePath(outputDirectory)


  start_time <- Sys.time()
  script_contents <- readr::read_file(system.file("scripts/StitchTiffChannels.py", package = "Flat2OMETiff"))
  script <- tempfile()

  #write the function definition to the temp file
  readr::write_file(script_contents, script)

    #convert the R vector to a JSON array to be read by python
  channelWhitelistJSON <- jsonlite::toJSON(channelWhitelist, auto_unbox = TRUE)
  print(channelWhitelistJSON)
  function_call <- paste0("StitchTiffChannels(channelsDirectory = '", tiffChannelDirectory,
                            "',outputDirectory = '", outputDirectory,
                            "',channelNames = '", channelWhitelistJSON,
                            "')")

    #write the function call with arguments to the end of the script and execute
  readr::write_file(function_call, script, append = TRUE)
  system2(reticulate::py_exe(), script)
  end_time <- Sys.time()
  print(paste0("Time to combine channels: ", end_time - start_time))

}
