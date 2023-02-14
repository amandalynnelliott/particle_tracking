# Amanda Elliott
# Institute for Modeling Plasmas and Cosmic Dust (IMPACT) - Boulder, CO
# Dust Lofted Tracking Script for Dust Lofted during CHARGING
# Nov and Dec 2022
# Last updated: 1/22/2023. 
# Based on code written for dust lofted during high voltage -- Version 20.
# Contact: amel3682@colorado.edu


# Process video/images from the Phantom V2512 high-speed camera.

# Detects particle movement and tracks particles between frames. 

# Data collected:
# --- (1) frame number
# --- (2) dust particle number
# --- (3) dust size 
# --- (4) velocity --> initial vertical velocity (v_iy = [v * sin(theta)] - only if v_y > 0)
# ------ v_x and v_y
# --- (5) surface potential phi(radius, velocity)

# Input parameters about video frames needed:
# --- (1) Frames per second that the video was recorded at
# --- (2) Pixels to microns. How many microns are in one pixel? -- For unit conversion.
# --- (3) If video is of dust sample in a perfect circle, need d1 (full diameter) and d2 (truncated diameter).
# ---------- Boolean if frames are perfectly horizontal or at an angle. If at an angle, need to calculate it (theta).
# ---------- theta = arccos(d2/d1)
 
# Imports -------------------------------------------------------------------
from PIL import Image, ImageDraw
import numpy as np
from math import floor, ceil, pi, sqrt, sin
import sys
import csv

import tracking_functions as tf


args = tf.command_line_args()

directory = args.directory
frames_start = args.start
frames_end = args.end

print(f"frames_start: {frames_start}")
print(f"frames_end: {frames_end}")
print(f"Number of frames: {frames_end - frames_start + 1}")

avg_brightness_sample_size = args.sample
movement_brightness_threshold = args.movement
dust_brightness_threshold = args.dust
t = args.circle


# Get values saved in video_text.txt file.
fps = 0
total_time = 0
A_focus = 0
lofting_rate_particle_counter = 0
camera_angle = 0 

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

    px_to_m = floats[1] * pow(10, -6)
    print(f"Pixels to meters: {px_to_m}")

    beam_voltage = floats[2]
    beam_current = floats[3]
    print(f"Charging electron beam at {beam_voltage} V and {beam_current} mA.")

    dust_density = floats[4]
    print(f"Dust density: {dust_density} %.")

    camera_angle = floats[5]  

    d2 = floats[6]
    d1 = floats[7]
    depth_of_field = floats[8]  

    theta = tf.calc_camera_angle_in_rad(camera_angle, d1, d2, depth_of_field)

    total_time = (frames_end - frames_start + 1) / fps
    print(f"Total time for frames: {total_time} seconds")

    f.close()
    

dt = 1/fps # Time in [sec] that elapses between frames.
density = tf.density
g = tf.g
epsilon_0 = tf.epsilon_0
vector_accent_color = tf.vector_accent_color
box_line_color = tf.box_line_color


# -------------------------------------------------------------------- 

def theo_phi(r, v, gamma = 1):
    '''
    Function: For calculating surface potential (phi) of dust particle.
    Input radius size, initial launch velocity, and gamma (= either 1 or 5).
    Reference: Anthony Carroll (2020), "Laboratory measurements of initial launch velocities of electrostatically 
            lofted dust on airless planetary bodies" 
    '''
    return (v * r / gamma) * np.sqrt(density / (3 * epsilon_0)) 

# -------------------------------------------------------------------
# -------- MAIN ----------------------------------------------------
# -------------------------------------------------------------------

width, height = tf.get_frame_dimensions(directory, frames_start)


# Brightness Maps -----------------------------------------------------

print()

brightness_map_3D = tf.load_brightness_maps(directory, frames_start, frames_end, width, height)

# -------- Create folders -----------------------------------------------------

tf.create_output_folders(directory)

#----------------------------------------------------------------------

# Generate a .csv file to write data and parameters to

print(" ")
column_spacing = "{:<8} {:<10} {:<15} {:<10} {:<16} {:<20} {:<15} {:<15} {:<20} {:<20} {:<15} {:<15}"

parameters_header = ['Directory', 'Frames per second (recorded)', 'Pixels to microns', 'Movement THRES','Dust THRES', 'Sample size (nxn)']
parameters = [directory, f'{fps}', f'{px_to_m / pow(10,-6)}', f'{movement_brightness_threshold}', f'{dust_brightness_threshold}', f'{avg_brightness_sample_size}']


parameters_2_header = ['Beam Voltage (V)', 'Beam Current (V)', 'Dust Density (%)' , 'Camera Angle (rad)']
parameters_2 = [f'{beam_voltage}', f'{beam_current}', f'{dust_density}', f'{theta}']

print("Checking for movement...") 
header = ['Frames', 'Particle', 'Location', 'Size [px]', 'Size [diam, um]', 'Velocity [cm/sec]', 'dv [cm/sec]', 'a [m/sec^2]', 'V_x [cm/sec]', 'V_y [cm/sec]',  'Height [px]', 'Phi [V]']
print(column_spacing.format(header[0], header[1], header[2], header[3], header[4], header[5], header[6], header[7], header[8], header[9], header[10], header[11]))

with open(directory + "/data.csv", "w") as c:
    writer = csv.writer(c)
    writer.writerow(parameters_header)
    writer.writerow(parameters)
    writer.writerow(parameters_2_header)
    writer.writerow(parameters_2)
    writer.writerow([''])
    writer.writerow(header)

# ---------------------------------------------------------------------------------------
# Check for movement via brightness changes.

lofted_particles = 0

brightness_difference = 0
avg_brightness_prev = 0
avg_brightness_current = 0 

particle_centers = []
sizes = []
velocities = []

particle_num = 0

for N in range(0, frames_end - frames_start + 1):  
    particles = []  # Individual particles contained within:  [[(),(),()],[(),()]] <-- 
    dv_s = []

    if (N==0): 
        particle_centers.append([])
        sizes.append([])

    else: 
        # Loop over image.
        for i in range(floor(avg_brightness_sample_size/2),width, avg_brightness_sample_size): 
            for j in range(floor(avg_brightness_sample_size/2),height, avg_brightness_sample_size):

                if(i < floor(avg_brightness_sample_size/2)): 
                    continue 
                if(i > width - ceil(avg_brightness_sample_size/2)):
                    continue
                if(j < floor(avg_brightness_sample_size/2)): 
                    continue 
                if(j > height - ceil(avg_brightness_sample_size/2)):
                    continue

                # Main algorithm for checking for movement ----------------------
                if(avg_brightness_sample_size == 1):
                    brightness_difference = brightness_map_3D[N][i][j] - brightness_map_3D[N-1][i][j]
                else: 

                    for x in range(i-floor(avg_brightness_sample_size/2), i+ceil(avg_brightness_sample_size/2)):
                        for y in range(j-floor(avg_brightness_sample_size/2), j+ceil(avg_brightness_sample_size/2)):
                            avg_brightness_prev += brightness_map_3D[N-1][x][y]


                    avg_brightness_prev /= avg_brightness_sample_size**2

                    for x in range(i-floor(avg_brightness_sample_size/2), i+ceil(avg_brightness_sample_size/2)):
                        for y in range(j-floor(avg_brightness_sample_size/2), j+ceil(avg_brightness_sample_size/2)):
                            avg_brightness_current += brightness_map_3D[N][x][y]


                    avg_brightness_current /= avg_brightness_sample_size**2

                brightness_difference = avg_brightness_current - avg_brightness_prev
                # -------------------------------------------------------------------

                if (brightness_difference > movement_brightness_threshold):

                    new_particle = True

                    for particle in particles: 

                        if(new_particle == False):
                            break
                        
                        for pixel in particle:
                            (p_x, p_y) = pixel

                            if((abs(i-p_x) <= 10) and (abs(j-p_y) <= 10)):
                                particle.append((i,j)) 
                                new_particle = False
                                break

                    if(new_particle == True):
                        particles.append([(i,j)])

        # Find the true center for each particle
        temp_centers = []
        temp_sizes = []

        for particle in particles:  # For finding the true center and also the size of each particle. 
            mid_x = 0
            mid_y = 0
    
            for pixel in particle:
                x,y = pixel
                mid_x += x
                mid_y += y

            mid_x /= len(particle)
            mid_y /= len(particle)

            temp_centers.append((mid_x,mid_y))

            # Flood-Fill Algorithm
            it_mid_x = int(mid_x)   # Iterable
            it_mid_y = int(mid_y)
            prevB = brightness_map_3D[N][it_mid_x][it_mid_y]
            
            size = tf.getSize(brightness_map_3D[N], width, height, it_mid_x, it_mid_y, prevB, dust_brightness_threshold)

            temp_sizes.append(size)

        
        particle_centers.append(temp_centers)
        sizes.append(temp_sizes)

        # Save images 
        if(len(particles) > 0):
            
            # Verify particle locations 
            if(N==1):
                first_frame = directory + "/frame" + str(frames_start) + ".jpg"
                img_first = Image.open(first_frame) 

                img_first.save(directory + "/lofted/frame" + str(frames_start) + ".jpg")

            filename = directory + "/frame" + str(N + frames_start) + ".jpg"
            img_N = Image.open(filename)
            img_N_copy = img_N.copy()
            draw_N = ImageDraw.Draw(img_N_copy)
            
            pairs = []

            for particle_center in particle_centers[N]: 
                (x,y) = particle_center
                draw_N.ellipse([x-t,y-t,x+t, y+t], fill=None, outline="red") #perfect circle

                smallest_dist = sys.maxsize
                smallest_index = -1

                # Find which particles between frames are closest together. 
                for i in range(len(particle_centers[N - 1])):
                    (xr,yr) = particle_centers[N-1][i]
                    dist = sqrt((x-xr)**2 + (y-yr)**2)
                
                    if(dist < smallest_dist):
                        smallest_dist = dist       # in units: pixels
                        smallest_index = i

                if(smallest_index != -1):
                    # Original velocity calculation from HV script.
                    # velocity = smallest_dist * fps * sin(theta) # in units: pixels/sec
                    pairs.append({"p1": particle_center, "p2": particle_centers[N-1][smallest_index], "dist": smallest_dist})
                    
            pairs_copy = list(pairs) # Will be removing items from pairs.

            for pair in pairs_copy:

                for other_pair in pairs_copy:

                    if(pair == other_pair): continue
                    elif((pair["p1"] == other_pair["p1"]) or (pair["p2"] == other_pair["p2"]) or (pair["p1"] == other_pair["p2"]) or (pair["p2"] == pair["p1"])):

                        if(pair["dist"] < other_pair["dist"]): 
                            if other_pair in pairs:
                                pairs.remove(other_pair)
                            break

            pairs_copy = list(pairs)

            for pair in pairs:
                draw_N.line(((pair["p2"][0], pair["p1"][1]),pair["p1"]), width=1, fill=box_line_color) # white bounding box
                draw_N.line(((pair["p1"][0], pair["p2"][1]),pair["p1"]), width=1, fill=box_line_color) # white bounding box
                tf.arrowedLine(draw_N, pair["p2"], pair["p1"], 1) # particle vector
                tf.arrowedLine(draw_N, pair["p2"], (pair["p2"][0], pair["p1"][1]), 1, vector_accent_color) # y-vector
                tf.arrowedLine(draw_N, pair["p2"], (pair["p1"][0], pair["p2"][1]), 1, vector_accent_color) # x-vector

                dy = (pair["p2"][1] - pair["p1"][1])              # height difference 
                v_x = (pair["p1"][0] - pair["p2"][0]) * fps       # x-velocity in px/sec
                v_y = (dy) * fps * sin(theta)    # y-velocity in px/sec. Corrected with camera angle.
                velocity = sqrt(v_x**2 + v_y**2) # total velocity in px/sec

                velocities.append({"position": pair['p1'], "position_prev": pair['p2'], "v_x" : v_x, "v_y" : v_y, "dy" : dy})
                
            img_N_copy.save(directory + "/lofted/frame" + str(N + frames_start) + ".jpg")

        # Get change in velocity, dv = v_current - v_prev
        for v_current in velocities:
            for v_prev in velocities:
                if(v_current["position_prev"] == v_prev["position"]):
                    
                    height_y = v_current["dy"] + v_prev["dy"]  # v_current = v_prev + dy 
                    v_prev["dy"] = 0    # Need to do this, otherwise dy values will be counted multiple times. 
                    v_current["dy"] = height_y
                    dv_s.append({"position": v_current["position"], "total_height": height_y, 
                            "v_y_current": v_current["v_y"], "v_y_prev": v_prev["v_y"],
                            "v_x_current": v_current["v_x"], "v_x_prev": v_prev["v_x"]})
                    

        lofting_rate_particle_counter += len(dv_s)

        # Print statments and save data to .csv file.
        particle_num = 0
        for particle in particles:
            particle_num += 1

            size_m = sizes[N][particle_num-1] * (px_to_m**2)  # Convert pixels to m^2
            size_r = sqrt(size_m/pi)
            volume = (4/3) * pi * pow(size_r, 3)  # particle volume in m^3
            size_diam = 2*size_r
            mass = density * volume     # particle mass in kg

            temp_velocity_x = 0
            temp_velocity_x_m = 0

            temp_velocity_m = 0

            temp_v_x_m = 0
            temp_v_y_m = 0
            temp_v_y_m_updated = 0

            temp_dv_m = 0

            lofted_theta = 0
            phi = 0

            H_lofted_current = 0
            H_lofted_prev = 0

            acceleration = 0  # acceleration 

            for delta_v in dv_s:
                if (delta_v["position"] == particle_centers[N][particle_num-1]):

                    # Calculate new v_y based on Anthony Carroll (2020): v_y_updated = sqrt(v_y**2 + (2*g*H_lofted))

                    # ----------- Current velocity -------------
                    H_lofted_current = delta_v["total_height"] #units: pixels

                    temp_v_x_current = delta_v["v_x_current"]
                    temp_v_x_m = temp_v_x_current * px_to_m    #units: m/s

                    temp_velocity_current_y_m = delta_v["v_y_current"] * px_to_m #units: m/sec
                    temp_v_y_m_updated = np.sqrt(((temp_velocity_current_y_m) ** 2) + (2 * g * abs(H_lofted_current * px_to_m)))

                    temp_velocity_m = np.sqrt(temp_v_x_m ** 2 + temp_v_y_m_updated ** 2) # Current launch velocity (m/sec)

                    # ----------- Previous velocity -------------
                    temp_v_x_prev_m = delta_v["v_x_prev"] * px_to_m         #units: m/sec
                    temp_velocity_prev_y_m = delta_v["v_y_prev"] * px_to_m  #units: m/sec

                    for delta_v_prev in dv_s:   # Get previous lofted height.
                        for particle_pos in particle_centers[N-1]:
                            if (delta_v_prev["position"] == particle_pos):
                                H_lofted_prev = delta_v_prev["total_height"] #units: pixels

                    temp_v_y_prev_updated_m = np.sqrt((temp_velocity_prev_y_m ** 2) + (2 * g * abs(H_lofted_prev * px_to_m)))
                    temp_velocity_prev_m = np.sqrt(temp_v_x_prev_m ** 2 + temp_v_y_prev_updated_m ** 2) # Previous launch velocity (m/sec)
                    
                    # ----------- Calculate acceleration -------------
                    temp_dv_m = temp_velocity_m - temp_velocity_prev_m
                    acceleration = temp_dv_m / dt #units: m/sec^2

                    # ----------- Calculate phi ----------------------
                    phi = theo_phi(size_r, temp_velocity_m)

            print(column_spacing.format(
                    "---" if(particle_num > 1) else f"{N - 1 + frames_start}-{N + frames_start}", 
                    "{:.0f}".format(particle_num), 
                    "({:.0f},{:.0f})".format(particle_centers[N][particle_num-1][0], particle_centers[N][particle_num-1][1]),
                    "{:.1f}".format(sizes[N][particle_num-1]),  
                    "{:.1f}".format(size_diam * pow(10,6)), 
                    "{:.1f}".format(temp_velocity_m * pow(10,2)),
                    "{:.1f}".format(temp_dv_m * pow(10,2)),
                    "{:.6f}".format(acceleration),
                    "{:.1f}".format(temp_v_x_m * pow(10,2)),
                    "{:.1f}".format(temp_v_y_m_updated * pow(10,2)),
                    "{:.4f}".format(H_lofted_current),
                    "{:.1f}".format(phi)
            ))

            data_row = ['---' if(particle_num > 1) else f"'{N-1 + frames_start}-{N + frames_start}",
                        particle_num,
                        "(" + str(particle_centers[N][particle_num-1][0]) + "," + str(particle_centers[N][particle_num-1][1]) + ")",
                        sizes[N][particle_num-1],
                        size_diam * pow(10,6),
                        temp_velocity_m * pow(10,2),
                        temp_dv_m * pow(10,2),
                        acceleration,
                        temp_v_x_m * pow(10,2),
                        temp_v_y_m_updated * pow(10,2),
                        H_lofted_current,
                        phi]

            with open(directory + "/data.csv", "a") as c:
                writer = csv.writer(c)
                writer.writerow(data_row)

        print(column_spacing.format(
                "--------",
                "---------",
                "--------------",
                "---------",
                "-----------",
                "---------------",
                "-------------",
                "-------------",
                "--------------",
                "--------------",
                "---------",
                "--------------",
                "--------------",
                "--------------",
                "---------"))

# Close .csv file that the data was written to. 
c.close()

# ---------------- USEFUL PLOT: Brightness as a 2D heat map -----------------

size_plt = args.sizes #read cmd param
if size_plt: tf.plot_brightness(size_plt, width, height, brightness_map_3D, frames_start, directory)

# --------------------------------------------------------------------------

print()

print("Done.")