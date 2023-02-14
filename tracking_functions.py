import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import math
from PIL import Image
import argparse
import sys
import os
import uproot 

# ------------------------------------------------------------------------
# Constants

density = 2.20 * pow(10,3) # density for SiO2 (silica) in units [kg/m^3]
g = 9.8 # gravitational const accel. m/s^2
epsilon_0 = 8.854 * 10**(-12)
vector_accent_color = (0,0,255)
box_line_color = (255,255,255)

d = 1.1 * pow(10, -2) # distance in meters between mesh and dust = 1.1 cm


# ------------------------------------------------------------------------

def command_line_args():
    '''
    Parse command line arguments.
    '''

    # Initialize parser
    parser = argparse.ArgumentParser()

    # Adding command line arguments
    parser.add_argument("directory", type=str,
                        help = "Directory that contains the frames as jpegs \
                                to be processed.")
    parser.add_argument("--start", type=int, required=True,
                        help = "Starting frame number.")
    parser.add_argument("--end", type=int, required=True,
                        help = "Ending frame number.")
    parser.add_argument("--movement", type=int, default = 30,
                        help = "Brightness threshold for movement detection.")
    parser.add_argument("--dust", type = int, default = 5,
                        help = "Brightness threshold to get particle size.")
    parser.add_argument("--sizes", type = str, 
                        help = "Input the range of frames 'start-end' to generate \
                                brightness maps for to check the flood fill algorithm \
                                size calculation. ")
    parser.add_argument("--sample", type = int, default = 2,
                        help = "Set a value >1 that will be the length of the \
                                side of the brightness averaging box.")
    parser.add_argument("--circle", type=int, default = 10,
                        help = "Radius in pixels of the red circle drawn \
                                to track the particle.")

    print('Command Line Arguments')
    print(sys.argv)

    return parser.parse_args()

def calc_camera_angle_in_rad(camera_angle_deg, d1, d2, depth_of_field):
    '''
    Calculate the angle of the camera in radians.
    If the angle is not 90 degrees, also print the lofting area that is in focus.
    '''
    if(camera_angle_deg == -1):
        theta = math.acos(d2/d1)   
    else: theta = camera_angle_deg * math.pi/180

    print(f"Camera is at angle: {math.degrees(theta)} degrees; {theta} radians.")

    if(camera_angle_deg != 90) and (d1 != 0):
        # Calculate in-focus area of dust sample if deposited in perfect circle..
        r = 1.64/2  # Diameter for each dust sample is 1.64 cm
        px_to_cm = 1.64 / d1
        depth_of_field_half = 0.5 * depth_of_field * px_to_cm * (depth_of_field/d1) # Convert from px to cm.
        h = r - depth_of_field_half 
        A_circle = math.pi * r**2
        A_segment = r**2 * math.acos((r-h)/r) - (r-h) * math.sqrt(2 * r * (h - h**2)) # https://www.omnicalculator.com/math/segment-area
        A_focus = A_circle - 2 * A_segment

        print(f"dust samples radius: {r} cm")
        # print(f"h: {h} cm")
        # print(f"A_circle: {A_circle} cm^2")
        # print(f"A_segment: {A_segment} cm^2")
        print(f"A_focus: {A_focus} cm^2")  # The area that is in focus. 

    return theta

def arrowedLine(draw, ptA, ptB, width=1, color=(0,255,0)):
    """
    Draw line from ptA to ptB with arrowhead at ptB.
    Function source: https://stackoverflow.com/questions/63671018/how-can-i-draw-an-arrow-using-pil
    """
    # Draw the line without arrows
    draw.line((ptA,ptB), width=width, fill=color)

    # Now work out the arrowhead
    # = it will be a triangle with one vertex at ptB
    # - it will start at 80% of the length of the line
    # - it will extend 8 pixels either side of the line
    x0, y0 = ptA
    x1, y1 = ptB
    # Now we can work out the x,y coordinates of the bottom of the arrowhead triangle
    xb = 0.80*(x1-x0)+x0
    yb = 0.80*(y1-y0)+y0

    # Work out the other two vertices of the triangle
    # Check if line is vertical
    if x0==x1:
       vtx0 = (xb-5, yb)
       vtx1 = (xb+5, yb)
    # Check if line is horizontal
    elif y0==y1:
       vtx0 = (xb, yb+5)
       vtx1 = (xb, yb-5)
    else:
       alpha = math.atan2(y1-y0,x1-x0)-90*math.pi/180
       a = 3*math.cos(alpha)
       b = 3*math.sin(alpha)
       vtx0 = (xb+a, yb+b)
       vtx1 = (xb-a, yb-b)

    # Now draw the arrowhead triangle
    draw.polygon([vtx0, vtx1, ptB], fill=color)


# -----------------------------------------------------------------------------
# Functions: For Flood Fill algorithm:

def getSize(screen, m, n, x, y, prevB, threshold):
    '''
    Returns the size of the dust particle in pixels. Uses a Flood Fill algorithm.
    Source: https://www.geeksforgeeks.org/flood-fill-algorithm/
    screen: brightness_map for that frame.
    (m, n): image size.
    (x,y): dust particle center.
    prevB: the brightness at (x,y).
    newB: value to replace brightness (termination condition).
    threshold: dust brightness. 
    '''
    # print("Getting size... ")

    def isValid(screen, m, n, x, y, prevB, newB, threshold):
        if (x < 0) or (x >= m) or (y < 0) or (y >= n) or (abs(screen[x][y] - prevB) > threshold) or (screen[x][y] == newB):
            return False
        return True

    newB = 500  # Termination condition. 

    queue = []

    # Append the position of starting pixel of the component
    queue.append([x, y])

    # Change the pixel with the new brightness
    screen[x][y] = newB
    
    size_count = 1

    # While the queue is not empty i.e. the whole component having prevB brightness
    # is not containing the newB brightness
    while queue:

        # Dequeue the front node
        currPixel = queue.pop()

        posX = currPixel[0]
        posY = currPixel[1]

        # Check if the adjacent pixels are valid
        if isValid(screen, m, n, posX + 1, posY, prevB, newB, threshold):

            # Color with newC if valid and enqueue
            screen[posX + 1][posY] = newB
            queue.append([posX + 1, posY])
            size_count += 1

        if isValid(screen, m, n, posX-1, posY, prevB, newB, threshold):
            screen[posX-1][posY] = newB
            queue.append([posX-1, posY])
            size_count += 1

        if isValid(screen, m, n, posX, posY + 1, prevB, newB, threshold):
            screen[posX][posY + 1] = newB
            queue.append([posX, posY + 1])
            size_count += 1

        if isValid(screen, m, n, posX, posY-1, prevB, newB, threshold):
            screen[posX][posY-1] = newB
            queue.append([posX, posY-1])
            size_count += 1

    return size_count

# -----------------------------------------------------------------------------

def get_frame_dimensions(directory, frames_start):
    '''
    Function: Get the dimensions of the first frame in the video.
    '''
    first_frame = directory + "/frame" + str(frames_start) + ".jpg" 
    with Image.open(first_frame) as image:
        width, height = image.size
        
    print("Image size: " + str(image.size))
    
    return width, height


def load_brightness_maps(directory, frames_start, frames_end, width, height):
    with uproot.open(directory + "/brightness_maps.root") as f:
        brightness_map_3D = f["brightness_maps"]["brightness_map"].array(library="np")
        brightness_map_3D = np.transpose(brightness_map_3D, (0,1,2))
        print("Brightness maps loaded.")
        return brightness_map_3D


def create_output_folders(directory):
    '''
    Created lofted and sizes folders in the directory. 
    The lofted folder will contain the output tracking frames.
    The sizes folder is for manually saving size results.
    '''
    parent_dir = directory
    path = os.path.join(parent_dir, "lofted")
    if (os.path.exists(path) == False):
        print("\'lofted\' folder created.")
        os.mkdir(path)
    else: print("\'lofted\' folder exists.")

    path = os.path.join(parent_dir, "lofted/sizes")
    if (os.path.exists(path) == False):
        print("\'sizes\' folder created.")
        os.mkdir(path)
    else: print("\'sizes\' folder exists.")


def plot_brightness(size_plt, width, height, brightness_map_3D, frames_start, directory):
    '''
    USEFUL PLOT: Brightness as a 2D heat map.
    Will plot both the original frame and the brightness map of that frame after flood fill.
    This is important for verifying that the size of the particles is being calculated correctly.
    '''
    print("Generating brightness maps to verify size calculations...")
    arr = size_plt.split('-') # ie: ['5', '8']
    start = int(arr[0])
    end = int(arr[1])

    xx, yy = np.meshgrid(range(width), range(height))
    for plot in range(start, end + 1):

        show = np.flipud(brightness_map_3D[plot - frames_start].transpose()) # orientate brightness map correctly

        fig, axs = plt.subplots(2, 1, sharex=False, sharey=False)
        
        plt.title(f'Frame {plot}')
        axs[0].set_aspect('equal')
        axs[0].pcolor(xx,yy, show, vmin=0, vmax=500)
        # plt.title(f'Particle in Frame {plot}')

        img = mpimg.imread(directory + '/frame' + str(plot) + '.jpg')
        axs[1].imshow(img)
        plt.show()

