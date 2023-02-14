import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import pprint
import scipy
import scipy.linalg
import argparse
import sys
import os

# Input a .csv file with the first column containing the lofted particle's locations in pixels.
# IMPORTANT: If there is a missing frame, put an "M" in that location in the .csv file.
# Have the tracked frame images in a "lofted" folder.
# Output is in the "angle" folder.
# Will print out lofted angle and initial velocity in px/sec.

parser = argparse.ArgumentParser()

# Adding command line arguments
parser.add_argument("directory", type=str,
                    help = "Directory that contains the frames as jpegs \
                            to be processed.")
parser.add_argument("--start", type=int, required=True,
                    help = "Starting frame number.")
parser.add_argument("--end", type=int, required=True,
                    help = "Ending frame number.")

args = parser.parse_args()

print(sys.argv)

directory = args.directory
start_frame = args.start
end_frame = args.end

with open(directory + "\\..\\" + "video_info.txt") as f:

    data = f.read().split()
    floats = []
    for elem in data:
        try:
            floats.append(float(elem))
        except ValueError:
            pass

    fps = floats[0]
    print(f"Frames per second (recorded): {fps}")

    dt = (1/fps)

    px_to_m = floats[1] * pow(10, -6)
    print(f"Pixels to meters: {px_to_m}")

    f.close()

print(f"frames_start: {start_frame}")
print(f"frames_end: {start_frame}")


# Create folder for results.
parent_dir = directory + "/lofted"
path = os.path.join(parent_dir, "angle")
if (os.path.exists(path) == False):
    print("\'angle\' folder created.")
    os.mkdir(path)
else: print("\'angle\' folder exists.")

t = 20  # Radius size for drawn circles.

skip = [0]
b = pd.read_csv(directory + '/locations.csv', header=None, skiprows=skip, skip_blank_lines=True)
locations = b[0]
print(locations)

X_vec = [0]*len(locations)
Y_vec = [0]*len(locations)
t_vec = [0]*len(locations)

skip_frames = []
missing_index = 0

for i in range(0,len(locations)):

    if locations[i] == "M":
        skip_frames.append(i + start_frame)
        missing_index = missing_index + 1
    else: 
        # Remove parentheses
        locations[i] = locations[i].replace('(', '')
        locations[i] = locations[i].replace(')', '')
    
        locations[i] = locations[i].split(',')
        X_vec[i - missing_index] = float(locations[i][0])
        Y_vec[i - missing_index] = float(locations[i][1])
        t_vec[i - missing_index] = i * dt

    print(locations[i])

for i in range(len(skip_frames)):
    X_vec.pop()
    Y_vec.pop()
    t_vec.pop()


print(f'X_vec: {X_vec}')
print(f'Y_vec: {Y_vec}')
print(f'skip_frames: {skip_frames}')
print(f't_fec: {t_vec}')


# Draw locations on frames to see parabola
t=3
for N in range(start_frame, end_frame + 1):

    if N not in skip_frames:

        filename = directory + "/lofted/frame" + str(N) + ".jpg"
        img_N = Image.open(filename)
        img_N_copy = img_N.copy()
        draw_N = ImageDraw.Draw(img_N_copy)

        skip_temp = 0

        for M in skip_frames:
            if N > M:
                skip_temp = skip_temp + 1

        for i in range(N - start_frame - skip_temp):

            x = X_vec[i]
            y = Y_vec[i]

            draw_N.ellipse([x-t,y-t,x+t, y+t], fill=None, outline="green")
            
        img_N_copy.save(directory + "/lofted/angle/frame" + str(N) + ".jpg")
print('Locations drawn.')

# Form matrix A for QR Factorization
# Of the form: A = [1, t_1, t_1^2, ... , 1, t_n, t_n^2]
temp = t_vec
power_value = 2

last_column = [pow(item, power_value) for item in temp]
A_temp = np.transpose([np.ones(len(t_vec)), t_vec, last_column])
print(A_temp)
# [ [], [], [] ] --> []
print("A better form:")
for m in A_temp:
    print(m)

# QR Factorization using Scipy library
A = np.array(A_temp)  # From the Wikipedia Article on QR Decomposition
Q, R = scipy.linalg.qr(A_temp, mode='economic')

print("A:")
pprint.pprint(A)

print("Q:")
pprint.pprint(Q)

print("R:")
pprint.pprint(R)

# Check that A = QR
# test_A = np.dot(Q,R)
# print("A = QR:")
# pprint.pprint(test_A)

# Find solution for x(t) = a + b * t + c * t^2
# x_soln = [a, b, c]
b_x = X_vec
x_soln = np.dot(scipy.linalg.inv(R), np.transpose(Q))
x_soln = np.dot(x_soln, b_x)
print(f'x_soln: {x_soln}')

# Find solution for y(t) = d + e * t + f * t^2
# y_soln = [d, e, f]
b_y = Y_vec
y_soln = np.dot(scipy.linalg.inv(R), np.transpose(Q))
y_soln = np.dot(y_soln, b_y)
print(f'y_soln: {y_soln}')

# Draw (x(t), y(t)) on last frame

filename = directory + "/lofted/angle/frame" + str(end_frame) + ".jpg"
img_N = Image.open(filename)
img_N_copy = img_N.copy()
draw_N = ImageDraw.Draw(img_N_copy)

t = 1
finer_dt = dt/10

for time in t_vec:
    for i in range(0,10,1):

        temp_t = time + (finer_dt * i)

        x = x_soln[0] + x_soln[1]*temp_t + x_soln[2] * temp_t**2
        y = y_soln[0] + y_soln[1]*temp_t + y_soln[2] * temp_t**2

        draw_N.ellipse([x-t,y-t,x+t, y+t], fill=None, outline="yellow")
    
img_N_copy.save(directory + "/lofted/angle/frame" + str(end_frame) + "_QRresult.jpg")
print('QR result drawn.')

# Lofted Angle Calculation
# Take derivative
init_time = 0
v_x = x_soln[1] + 2 * x_soln[2] * init_time
v_y = y_soln[1] + 2 * y_soln[2] * init_time

theta = np.arctan2(v_y,v_x)
theta_deg = theta * 180 / np.pi
print(f'theta (rad): {theta}')
print(f'theta (deg): {theta_deg}')

initial_velocity = np.sqrt(v_x **2 + v_y **2)
print(f"initial velocity (px/sec): {initial_velocity}")

initial_velocity_cm = initial_velocity * px_to_m * 100
print(f"initial velocity (cm/sec): {initial_velocity_cm}")

initial_velocity_cm_x = v_x * px_to_m * 100
print(f"initial v_x (cm/sec): {initial_velocity_cm_x}")

initial_velocity_cm_y = v_y * px_to_m * 100
print(f"initial v_y (cm/sec): {initial_velocity_cm_y}")

# Check Angle
r = 100
filename = directory + "/lofted/angle/frame" + str(end_frame) + "_QRresult.jpg"
img_N = Image.open(filename)
img_N_copy = img_N.copy()
draw_N = ImageDraw.Draw(img_N_copy)
ptA = (X_vec[0], Y_vec[0])
ptB = (X_vec[0] + r, Y_vec[0])
draw_N.line((ptA,ptB), width=2, fill='red')

ptB = (ptA[0] + r * np.cos(theta), ptA[1] + (r * np.sin(theta)))
draw_N.line((ptA,ptB), width=2, fill='red')

text = "{:.1f}".format((-1) * theta_deg)+ ' degrees'
draw_N.text((ptA[0] + r/2, ptA[1] - r/4), text, align ="left") 

img_N_copy.save(directory + "/lofted/angle/frame" + str(end_frame) + "_QRresult_angle.jpg")

# For checking with alternative method
total_x = X_vec[len(X_vec) - 1] - X_vec[0]
total_y = Y_vec[len(Y_vec) - 1] - Y_vec[0]
t_f = t_vec[len(t_vec) - 1]

print(f"x displacement (px): {total_x}")
print(f"y displacement (px): {total_y}")
print(f"final time (s): {t_f}")