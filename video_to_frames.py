# Amanda Elliott
# June 2022
# Updated 2/9/2023

# Run this script to capture the frames of a video file as jpegs.


# IMPORTS ---------------------------------------------------------------
import sys 
import cv2
import numpy as np
from PIL import Image, ImageEnhance
import os
import argparse
import uproot

# COMMAND LINE ARGUMENTS ------------------------------------------------

# Initialize parser
parser = argparse.ArgumentParser()

# Adding command line arguments
parser.add_argument("directory", type=str,
                    help = "Directory that contains the frames as jpegs \
                            to be processed.")
parser.add_argument("video_name", type=str,
                    help = "Name of video file.")
parser.add_argument("--contrast", type = float, default = 1.0,
                    help = "Floating point. Factor 1.0 returns copy \
                        of original image. Lower # means less contrast. \
                        Higher # means more.")

args = parser.parse_args()
contrast = args.contrast
print(contrast)

print(sys.argv)

directory = args.directory
video_name = args.video_name

# -----------------------------------------------------------------------
# Function: Pass in the folder where the video file is located, and the name of
# the video file. The same folder will be filled with jpegs (individual frames)
# frpm the video. Returns the number of frames.  
def video_to_frames(directory, video_name):
    vidcap = cv2.VideoCapture(directory + "/" + video_name)
    success,image = vidcap.read()
    count = 0

    while success:
        cv2.imwrite(directory + "/" + "frame%d.jpg" % count, image)     # save frame as JPEG file      
        success,image = vidcap.read()
        # print('Read a new frame: ', success)
        count += 1

    return count


# ---------------------------------------------------------------------
# Function: Get video fps
def video_fps(directory, video_name):
    vidcap = cv2.VideoCapture(directory + "/" + video_name)
    fps = int(vidcap.get(cv2.CAP_PROP_FPS))

    return fps

# -----------------------------------------------------------------------
# Function: Convert image to brightness map array (2D)
def generate_brightness_arr(img, width, height):
    img = img.convert('RGB')
    brightness_map = np.zeros(image.size, dtype=float)

    for i in range(width):
        for j in range(height):
            pixelRGB = img.getpixel((i,j))
            R,G,B = pixelRGB
            pixel_brightness = sum([R,G,B])/3
            brightness_map[i][j] = pixel_brightness 

    return brightness_map


# MAIN -----------------------------------------------------------------------

# Import video file(s) and convert them to image frames

print("Processing video: ")
print("Getting frames...")

frames = video_to_frames(directory, video_name) 
print(f"Frames: {frames}")

fps = video_fps(directory, video_name)
print(f"Frames per second: {fps}")


# Get dimensions from first frame.
first_frame = directory + "/frame0.jpg" 
with Image.open(first_frame) as image:
    width, height = image.size

print("Image size: " + str(image.size))

if (contrast != 1):
    print("Increasing contrast...")

    # Create contrast folder
    path = os.path.join(directory, "original")
    if (os.path.exists(path) == False):
        print("\'original\' folder created.")
        os.mkdir(path)
    else: print("\'original\' folder exists.")

    for N in range(0, frames):

        filename = directory + "/frame" + str(N) + ".jpg"
        img = Image.open(filename)
        img.save(directory + "/original/frame" + str(N) + ".jpg")

        # delete original
        # os.remove(filename)
        img_contrast = ImageEnhance.Contrast(img)
        img_contrast.enhance(contrast).save(directory + "/frame" + str(N) + ".jpg")  


# Generate brightness maps as a 3D array. 
print("Generating brightness maps... ")
brightness_map_3D = []
for i in range(0, frames):
    filename = directory + "/frame" + str(i) + ".jpg"
    img = Image.open(filename)
    arr = generate_brightness_arr(img, width, height) 
    brightness_map_3D.append(arr)
    print('Brightness map for frame' + str(i) + " generated.")

# Save brightness maps to root file. 
print("Saving to root file...")
with uproot.recreate(directory + "/brightness_maps.root") as f:
    f["brightness_maps"] = {"brightness_map": np.array(brightness_map_3D)}
    f["image_specs"] = {"width": np.array([width]), "height": np.array([height]), "frames": np.array([frames])}
print('Saved.')
