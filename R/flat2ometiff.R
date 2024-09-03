#' @title Flat2Matrix
#'
#' @description
#' Temporary function to set up R package
#' @param echo boolean to control responses printed (or is it???)
#' @export

Flat2Matrix <- function(expressionCsvPath = NULL,
                        polygonCsvPath = NULL,
                        basePath = "./",
                        channelNamesAreExpressionCSVColumnNames = TRUE) {
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

  for (channel in channel_names) {
    #write blank image with correct global dimensions
    sparseMatrix <- .InitalizeSparseMatrix(x_min, x_max, y_min, y_max)

    #this sets up a vector of cells within FOVs that can be strsplit via:
    #FOV ID = strsplit(unique_cell_ids, split ="_")[2]
    #Cell ID = strsplit(unique_cell_ids, split ="_")[4]
    unique_cells <- paste0("fov_", flat_file_list$expression$fov, "_cellID_", flat_file_list$expression$cell_ID)

    #iterate over cells and fill the convex hulls (defined by the polygon flat files) with expression data
    for (cell in unique_cells[1:3]){
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

      #use the original coordinates directly in the sparse matrix
      i_indices <- cell_sub_image$x
      j_indices <- cell_sub_image$y

      #adjust indices to start from 1 within the bounding box
      i_adjusted <- i_indices - min(i_indices) + 1
      j_adjusted <- j_indices - min(j_indices) + 1

      #filter out the zero values to keep the matrix sparse
      non_zero_indices <- which(inside != 0)

      #define the sparse matrix just for the bounding box of the cell
      sub_region_nrows <- max(i_indices) - min(i_indices) + 1
      sub_region_ncols <- max(j_indices) - min(j_indices) + 1
      i_adjusted <- i_indices[non_zero_indices] - min(i_indices) + 1
      j_adjusted <- j_indices[non_zero_indices] - min(j_indices) + 1

      #calculate global indices for updating the larger matrix
      global_j_indices <- j_adjusted + min(j_indices)
      global_i_indices <- i_adjusted + min(i_indices)

      #check for and handle zero values
      #the zeroes in the sparse matrix will yield NA, since they're still sparse, so . != 0 yields NA
      valid_indices <- !is.na(global_i_indices) & !is.na(global_j_indices)


      #filter valid indices
      valid_global_i_indices <- global_i_indices[valid_indices]
      valid_global_j_indices <- global_j_indices[valid_indices]
      valid_inside_values <- inside[non_zero_indices][valid_indices]
      valid_i_adjusted <- i_adjusted[valid_indices]
      valid_j_adjusted <- j_adjusted[valid_indices]

      #create the sparse matrix
      cell_sparse_matrix <- Matrix::sparseMatrix(
        i = valid_i_adjusted,
        j = valid_j_adjusted,
        x = valid_inside_values,
        dims = c(sub_region_nrows, sub_region_ncols)
      )

      sparseMatrix[valid_global_i_indices, valid_global_j_indices] <-
        sparseMatrix[valid_global_i_indices, valid_global_j_indices] +
        cell_sparse_matrix[valid_i_adjusted, valid_j_adjusted]
    }
    if (channelNamesAreExpressionCSVColumnNames) {
      channel_name <- channel
    } else {
      channel_name <- which(channel_names == channel)
    }
    #write the channel's sparse matrix to disk
    Matrix::writeMM(sparseMatrix, file = paste0(basePath, "/sparse_matrix_", channel_name, ".mtx"))
  }
}

.InitalizeSparseMatrix <- function(x_min, x_max, y_min, y_max) {
  sparse_matrix <- Matrix::Matrix(0,
                 nrow = length(seq(x_min, x_max, 1)),
                 ncol = length(seq(y_min, y_max, 1)),
                 sparse = TRUE)
  return(sparse_matrix)
}


.ReadFlatFiles <- function(expressionCsvPath = NULL,
                          polygonCsvPath = NULL) {
  expression <- data.table::fread(expressionCsvPath)
  polygons <- data.table::fread(polygonCsvPath)
  return(list(expression = expression, polygons = polygons))
}

