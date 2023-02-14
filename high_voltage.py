# Amanda Elliott
# Institute for Modeling Plasmas and Cosmic Dust (IMPACT) - Boulder, CO
# Dust Lofted Tracking Script
# June and July 2022
# Last updated: 10/15/2022. Version 20.
# Contact: amel3682@colorado.edu


# Process video/images from the Phantom V2512 high-speed camera.

# Detects particle movement and tracks particles between frames. 

# Data collected:
# --- (1) frame number
# --- (2) dust particle number
# --- (3) dust size 
# --- (4) velocity --> initial vertical velocity (v_iy = [v * sin(theta)] - only if v_y > 0)
# --- (5) dust charge

# Input parameters about video frames needed:
# --- (1) Frames per second that the video was recorded at
# --- (2) Pixels to microns. How many microns are in one pixel? -- For unit conversion.
# --- (3) If video is of dust sample in a perfect circle, need d1 (full diameter) and d2 (truncated diameter).
# ---------- Boolean if frames are perfectly horizontal or at an angle. If at an angle, need to calculate it (theta).
# ---------- theta = arccos(d2/d1)
# --- (4) The voltage at which the particle lofted at, in kV. 
 
# Imports -------------------------------------------------------------------
from PIL import Image, ImageDraw
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


# Get the frames per second saved by video_to_frames.py
fps = 0
total_time = 0
A_focus = 0
lofting_rate_particle_counter = 0
camera_angle_not_90 = 0 
d = tf.d 

with open(directory + "/video_HV_lofted.txt") as f:
    data = f.read().split()
    floats = []
    for elem in data:
        try:
            floats.append(float(elem))
        except ValueError:
            pass
    
    V = floats[0]*1000 # Voltage [V] when particle lofting occured. 
    E = V/d # Electric field generated [V/m]
    print(f"Lofted at {V} V")

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
    

# Constants
dt = 1/fps # Time in [sec] that elapses between frames.
density = tf.density

# -------------------------------------------------------------------
# -------- MAIN ----------------------------------------------------
# -------------------------------------------------------------------

width, height = tf.get_frame_dimensions(directory, frames_start)

# Brightness Maps -----------------------------------------------------

print()

brightness_map_3D = tf.load_brightness_maps(directory, frames_start, frames_end, width, height)


# -------- Create folders -----------------------------------------------------

tf.create_output_folders(directory)

# Generate a .csv file to write data and parameters to

print(" ")
column_spacing = "{:<8} {:<10} {:<15} {:<10} {:<12} {:<16} {:<20} {:<20} {:<15} {:<15} {:<10}"

parameters_header = ['Directory', 'Frames per second (recorded)', 'Pixels to microns', 'Movement THRES','Dust THRES', 'Sample size (nxn)']
parameters = [directory, f'{fps}', f'{px_to_m * pow(10,-6)}', f'{movement_brightness_threshold}', f'{dust_brightness_threshold}', f'{avg_brightness_sample_size}']


parameters_2_header = ['Beam Voltage (V)', 'Beam Current (V)', 'Dust Density (%)' , 'Lofted Voltage (V)', 'Camera Angle (rad)']
parameters_2 = [f'{beam_voltage}', f'{beam_current}', f'{dust_density}',f'{V}', f'{theta}']

header = ['Frames', 'Particle', 'Location', 'Size [px]', 'Size [um^2]', 'Size [diam, um]', 'Velocity [px/sec]', 'Velocity [cm/sec]', 'dv [cm/sec]', 'a [m/sec^2]', 'q [C] 10^-14']
print(column_spacing.format(header[0], header[1], header[2], header[3], header[4], header[5], header[6], header[7], header[8], header[9], header[10]))

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
print("Checking for movement...") 
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
                    velocity = smallest_dist * fps * sin(theta) # in units: pixels/sec
                    pairs.append({"p1": particle_center, "p2": particle_centers[N-1][smallest_index], "dist": smallest_dist, "velocity": velocity})
                    
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
                velocities.append({"position": pair['p1'], "position_prev": pair['p2'], "velocity": pair['velocity']})
                tf.arrowedLine(draw_N, pair["p2"], pair["p1"], 1)


            # img_copy.show()
            img_N_copy.save(directory + "/lofted/frame" + str(N + frames_start) + ".jpg")

        # Get dv = v_current - v_prev
        for v_current in velocities:
            for v_prev in velocities:
                if(v_current["position_prev"] == v_prev["position"]):
                    dv = v_current["velocity"] - v_prev["velocity"]
                    dv_s.append({"position": v_current["position"], "dv": dv})

                    # print(f'lofting_rate_particle_current (dv): {lofting_rate_particle_counter}')
                    # lofting_rate_particle_counter -= 1
                    # break
        lofting_rate_particle_counter += len(dv_s)

        #Print statments
        particle_num = 0
        for particle in particles:
            particle_num += 1

            size_m = sizes[N][particle_num-1] * (px_to_m**2)  # Convert pixels to m^2
            size_r = sqrt(size_m/pi)
            volume = (4/3) * pi * pow(size_r, 3)  # particle volume in m^3
            size_diam = 2*size_r
            mass = density * volume     # particle mass in kg

            temp_velocity = 0
            temp_velocity_m = 0

            temp_dv = 0
            temp_dv_m = 0
            acceleration = 0  # acceleration 
            charge = 0

            for velocity in velocities:
                if (velocity["position"] == particle_centers[N][particle_num-1]):
                    temp_velocity = velocity["velocity"]
                    temp_velocity_m = temp_velocity * px_to_m

            for delta_v in dv_s:
                if (delta_v["position"] == particle_centers[N][particle_num-1]):
                    temp_dv = delta_v["dv"]
                    temp_dv_m = temp_dv * px_to_m #units: m/s
                    acceleration = temp_dv_m / dt #units: m/sec^2
                    g = 9.8 # gravitational const accel. m/s^2
                    charge = ((mass * acceleration + (mass * g))/ E)  # in Coulombs

            print(column_spacing.format(
                    "---" if(particle_num > 1) else f"{N - 1 + frames_start}-{N + frames_start}", 
                    "{:.0f}".format(particle_num), 
                    "({:.0f},{:.0f})".format(particle_centers[N][particle_num-1][0], particle_centers[N][particle_num-1][1]),
                    "{:.1f}".format(sizes[N][particle_num-1]), 
                    "{:.1f}".format(size_m * pow(10,12)), 
                    "{:.1f}".format(size_diam * pow(10,6)),
                    "{:.1f}".format(temp_velocity), 
                    "{:.1f}".format(temp_velocity_m * pow(10,2)),
                    "{:.1f}".format(temp_dv_m * pow(10,2)),
                    "{:.6f}".format(acceleration),
                    "{:.6f}".format(charge * pow(10,14))))

            data_row = ['---' if(particle_num > 1) else f"'{N-1 + frames_start}-{N + frames_start}",
                        particle_num,
                        "(" + str(particle_centers[N][particle_num-1][0]) + "," + str(particle_centers[N][particle_num-1][1]) + ")",
                        sizes[N][particle_num-1],
                        size_m * pow(10,12),
                        size_diam * pow(10,6),
                        temp_velocity,
                        temp_velocity_m * pow(10,2),
                        temp_dv_m * pow(10,2),
                        acceleration,
                        charge * pow(10,14)]

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
                "-------------------",
                "-------------------",
                "--------------",
                "--------------",
                "---------"))

c.close()

    # print("Max brightness difference: " + str(max_brightness_difference))

# ---------------- USEFUL PLOT: Brightness as a 2D heat map -----------------

size_plt = args.sizes #read cmd param
if size_plt: tf.plot_brightness(size_plt, width, height, brightness_map_3D, frames_start, directory)
    
# --------------------------------------------------------------------------

print("Done.")
