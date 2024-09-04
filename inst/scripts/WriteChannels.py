import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path
import tifffile, os

def WriteChannelTiff(sparseMatrixFile, outputDirectory):
  #check if the directory exists and create it if it doesn't
  os.makedirs(outputDirectory, exist_ok=True)
  
  #read the sparse matrix, get the channel's name, and convert to csr matrix
  channel_name = sparseMatrixFile.split("_", 2)
  channel_name = channel_name[2].split(".", 1)[0]
  sparseMatrix = mmread(sparseMatrixFile)
  sparseMatrix = csr_matrix(sparseMatrix)

  #get the dimensions of the matrix
  height, width = sparseMatrix.shape

  #create a new tiff file for the channel
  output_file = outputDirectory + "/" + channel_name + ".tiff"

  #open the tiff file for writing
  with tifffile.TiffWriter(output_file, bigtiff=True) as tiff:
    #write the matrix data into the image row by row
    for row in range(height):
        #read the row from the sparse matrix
        row_data = sparseMatrix[row].toarray().flatten()
        
        #clip the values to be within the range 0 to 255
        row_data = np.clip(row_data, 0, 255)
        
        # cnvert the row's data to uint8
        row_data = row_data.astype('uint8')
        
        #reshape the row data to match the image width
        row_image = row_data.reshape(1, width)
        
        #write the row to the TIFF file
        tiff.write(row_image, contiguous=True)
