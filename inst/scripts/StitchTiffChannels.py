import numpy as np
from pathlib import Path
import tifffile, os, json
import xml.etree.ElementTree as ET

def StitchTiffChannels(channelsDirectory, outputDirectory, channelNames):
    #ensure the output directory exists, and creat it if not
    os.makedirs(outputDirectory, exist_ok=True)
    
    #convert channelNames from the input JSON string to a python list
    channelNames = json.loads(channelNames)
    
    #define the output OME-TIFF file name based on the channel names
    merged_output_file = os.path.join(outputDirectory, "_".join(channelNames) + ".ome.tiff")
    
    #list to hold the channel images
    channel_images = []
    image_shape = None

    #loop through each channel and read the corresponding TIFF file
    for channel_name in channelNames:
        found = False
        for file in os.listdir(channelsDirectory):
            if channel_name in file and (file.endswith('.tiff') or file.endswith('.tif')):
                file_path = os.path.join(channelsDirectory, file)
                print(f"Processing channel '{channel_name}' from file: {file_path}")
                
                #read the single-channel TIFF image
                img = tifffile.imread(file_path)
                print(f"Image shape for channel '{channel_name}': {img.shape}")
                
                #ensure the image shape is consistent across channels
                if image_shape is None:
                    image_shape = img.shape
                elif img.shape != image_shape:
                    raise ValueError(f"Image shape mismatch: expected {image_shape}, but got {img.shape} for channel {channel_name}")
                
                #add this channel's image to the list
                channel_images.append(img)
                found = True
                break
        
        if not found:
            raise FileNotFoundError(f"Channel '{channel_name}' not found in directory.")

    #stack the channel images into a multi-channel image (C, Y, X)
    multi_channel_img = np.stack(channel_images, axis=0)  # Shape: (C, Y, X)
    print(f"Final multi-channel image shape: {multi_channel_img.shape}")
    
    #prepare metadata for OME-TIFF
    metadata = {
        'axes': 'CYX',  #3D with channel axis
        'Channel': [{'Name': channel} for channel in channelNames],
        'Pixels': {
            'SizeC': len(channelNames),  #number of channels
            'SizeY': multi_channel_img.shape[1],  #Y
            'SizeX': multi_channel_img.shape[2],  #X
            'Type': 'uint8' if multi_channel_img.dtype == np.uint8 else 'uint8'
        }
    }
    
    #Write the stacked multi-channel as OME with metadata (and compression)
    tifffile.imwrite(merged_output_file, multi_channel_img, photometric='minisblack', metadata=metadata, bigtiff=True, compression='zlib')

        #confirm the size of the written image and the number of stacks
    with tifffile.TiffFile(merged_output_file) as tif:
        #read metadata for sizes
        image_metadata = tif.pages[0].tags
        image_width = image_metadata['ImageWidth'].value
        image_length = image_metadata['ImageLength'].value
        print(f"Written image dimensions: {image_length} x {image_width}")
        
        #extract and parse the OME-XML metadata
        ome_xml = image_metadata['ImageDescription'].value
        root = ET.fromstring(ome_xml)
        
        #print the name of each channel to ensure the channels are written correctly. 
        print("Channel names:")
        for channel in root.findall(".//{http://www.openmicroscopy.org/Schemas/OME/2016-06}Channel"):
            print(channel.get('Name'))
        
        #count the number of stacks/channeles and report
        num_stacks = len(tif.pages)
        print(f"Number of stacks (pages) in the TIFF file: {num_stacks}")

    print(f"Multi-channel OME-TIFF saved as: {merged_output_file}")

