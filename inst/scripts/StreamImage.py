from PIL import Image
import numpy as np
from scipy.io import mmread

def StreamImage(sparseMatrixFile):
  sparseMatrix = mmread(sparseMatrixFile)

  #matrix dimensions
  image_height, image_width = sparseMatrixshape 
  #Create an image object
  image = Image.new('L', (image_height, image_width))

  #process in chunks (chunk dimensions = the length of the image x 1000)
  chunk_height = 1000 
  chunk_width = image_width
  
  #loop through the image in chunks
  for y in range(0, image_height, chunk_height):
    #determine the actual height of the chunk (may be smaller for the last chunk)
    current_chunk_height = min(chunk_height, image_height - y)
    
    #create an empty chunk
    chunk = np.zeros((current_chunk_height, chunk_width), dtype=np.uint8)
    
    #get the indices for the current chunk
    mask = (sparseMatrix.row >= y) & (sparseMatrix.row < y + current_chunk_height)
    chunk_rows = sparseMatrix.row[mask] - y
    chunk_cols = sparseMatrix.col[mask]
    chunk_data = sparseMatrix.data[mask]
    
    #place the data into the chunk
    chunk[chunk_rows, chunk_cols] = chunk_data
    
    #convert the chunk to an image
    chunk_image = Image.fromarray(chunk)

    #paste the chunk into the final image at the correct position
    image.paste(chunk_image, (0, y))

  #save the final image
  image.save("non_square_large_image.png")
