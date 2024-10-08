% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flat2ometiff.R
\name{WriteTiffChannels}
\alias{WriteTiffChannels}
\title{WriteTiffChannels}
\usage{
WriteTiffChannels(
  sparseMatrixDirectory = "./sparse_matrix_files/",
  outputDirectory = "./tiff_channels/"
)
}
\arguments{
\item{sparseMatrixDirectory}{Path to the directory containing the sparse matrices.}

\item{outputDirectory}{Path to the directory where the tiff files will be written.}
}
\description{
This function writes the sparse matrices to tiff files using the python tifffile library.
}
