import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
import tifffile, os

def WriteChannelTiff(sparseMatrixFile, outputDirectory):
    #check if the directory exists and create it if it doesn't
    os.makedirs(outputDirectory, exist_ok=True)
    
    #parse the file path (the channel name will be in the last split on _ characters)
    channel_name = sparseMatrixFile.split("_")[-1].split(".")[0]
    print(f"Channel name: {channel_name}")
    
    #read the sparse matrix from R and convert to csr_matrix (different kind of sparse matrix)
    sparseMatrix = mmread(sparseMatrixFile)
    sparseMatrix = csr_matrix(sparseMatrix)
    
    #convert the sparse matrix to a supported data type (e.g., uint8). uint16 will double the file size of the matrix.
    sparseMatrix = sparseMatrix.astype(np.uint8)
    
    #get the dimensions of the matrix and print size
    height, width = sparseMatrix.shape
    print(f"Original matrix dimensions: {height} x {width}")

    #create a memory-mapped file for the numpy array, since densifying the matrix requires like ~60GB of RAM
    memmap_file = os.path.join(outputDirectory, f"{channel_name}.npy")
    memmap_array = np.memmap(memmap_file, dtype=np.uint8, mode='w+', shape=(height, width))

    #write the sparse matrix data to the memory-mapped file
    memmap_array[:] = sparseMatrix.toarray()
    memmap_array.flush()

    #create a new TIFF file for the channel
    output_file = os.path.join(outputDirectory, f"{channel_name}.ome.tiff")

    #write the entire image from the memory-mapped file, using compression
    tifffile.imwrite(output_file, memmap_array, photometric='minisblack', bigtiff=True, compression ='zlib')

    print(f"Finished writing TIFF to {output_file}")

    #confirm the size of the written image and the number of stacks (should match the original matrix dimensions, and number of stacks should be 1)
    with tifffile.TiffFile(output_file) as tif:
        # Access the metadata
        image_metadata = tif.pages[0].tags
        image_width = image_metadata['ImageWidth'].value
        image_length = image_metadata['ImageLength'].value
        print(f"Written image dimensions: {image_length} x {image_width}")

        # Count the number of stacks (pages)
        num_stacks = len(tif.pages)
        print(f"Number of stacks (pages) in the TIFF file: {num_stacks}")
    #close the memmap file
    memmap_array._mmap.close()
    #remove the memory mapped file 
    os.remove(memmap_file)

