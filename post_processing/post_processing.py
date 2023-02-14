import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

epsilon_0 = 8.854 * 10**(-12)
density = 2.20 * pow(10,3) # density for SiO2 (silica) in units [kg/m^3]

# -----------------------------------------------------------------------------
# Functions: For curve-filling and theoretical calculations. 

def theo_velocity(r, phi): 
    return ((phi / r) * np.sqrt((3 * epsilon_0 / density)))

def theo_phi(r, v, gamma=1):
    return (v * r / gamma) * np.sqrt(density / (3 * epsilon_0)) 


# -----------------------------------------------------------------------------
# Functions: Types of plots.


# -----------------------------------------------------------------------------
# Function: Generate plots for charging script. 
def plotter(sizes_diam, vy, vx, velos, thetas, phis):

    print(" ")
    print("Plotting... ")

    velos = velos / 100
    vy = vy/100
    sizes_r = sizes_diam/2

    print(f'sizes_r: {sizes_r}')
    # QV = 0.0668337^2/2 = 0.00223 ug m^2/s^2 = 2.23 * 10^-12 J
    # fit curve (of type func()) to samp data
    popt, pcov = curve_fit(theo_velocity, sizes_r, velos)  # popt is phi in microVolts
    x_sample = np.linspace(0.1, 75, 250)
    print(f'popt{popt * 10 ** (-6)}')  # prints phi in Volts
    # plot scatter and curve fit


    fig = plt.figure()
    # ax1 = fig.add_subplot(111)

    # ax1.scatter(sizes_r, velos, s=10, c='b', marker="s", label='first')
    # ax1.scatter(Ksizes_r, Kvelos, s=10, c='r', marker="o", label='second')

    scatter = plt.scatter(sizes_r, velos, c='k', label='Amanda_data')
    # corrected = plt.plot(x_sample, theo_velocity(x_sample, 6 * 22 * 10**6), color = 'blue' ) #, label=r'$QV = 0.0022$')
    # gamma2 = plt.plot(x_sample, theo_velocity(x_sample, 4 * 22 * 10 ** 6), color='green')
    gamma2 = plt.plot(x_sample, theo_velocity(x_sample, 5 * 22 * 10 ** 6), color = 'blue')
    theory = plt.plot(x_sample, theo_velocity(x_sample, 1 * 22 * 10 ** 6), color = 'red')
    #theory_lower = plt.plot(x_sample, func(x_sample, np.sin(30. * np.pi / 180.) * 22 * 10**6 / np.sin(50. * np.pi / 180.)), 'r--', dashes = (5, 10))
    #theory_upper = plt.plot(x_sample, func(x_sample, np.sin(70. * np.pi / 180.) * 22 * 10**6 / np.sin(50. * np.pi / 180.)), 'r--', dashes = (5, 10))
    #plt.fill_between(x_sample, func(x_sample, np.sin(30. * np.pi / 180.) * 22 * 10**6 / np.sin(50. * np.pi / 180.)), func(x_sample, np.sin(70. * np.pi / 180.) * 22 * 10**6 / np.sin(50. * np.pi / 180.)), color ='red', alpha=0.2)
    #corrected_lower = plt.plot(x_sample, func(x_sample, 4 * np.sin(30. * np.pi / 180.) * 22 * 10**6 / np.sin(50. * np.pi / 180.)), 'b--', dashes =(5, 10))
    #corrected_upper = plt.plot(x_sample, func(x_sample, 4 * np.sin(70. * np.pi / 180.) * 22 * 10 ** 6 / np.sin(50. * np.pi / 180.)), 'b--', dashes =(5, 10))
    #plt.fill_between(x_sample, func(x_sample, 4 * np.sin(30. * np.pi / 180.) * 22 * 10 ** 6 / np.sin(50. * np.pi / 180.)), func(x_sample, 4 * np.sin(70. * np.pi / 180.) * 22 * 10 ** 6 / np.sin(50. * np.pi / 180.)),color='blue', alpha=0.2)

    sml = 12
    med = 14
    lrg = 16

    plt.rc('xtick', labelsize = med)
    plt.rc('ytick', labelsize = med)
    plt.rc('legend', fontsize = sml)

    plt.legend(['Amanda\'s Data', '$\gamma$ = 5', '$\gamma$ = 1'])
    plt.xlabel('Dust Radius ($\mu m$)', fontsize = 20)
    plt.xticks(fontsize = 16)
    plt.ylabel('Launch Velocity ($m/s$)', fontsize = 20)
    plt.yticks(fontsize=16)
    plt.grid(True)
    x = sizes_r
    y = velos

    xerr = 0.2          # Fix this!
    yerr = [y*0, y*0.2] # Fix this!
    # plt.errorbar(x, y, yerr, xerr, fmt='none', ecolor='black', elinewidth=0.5)
    
    xtheory = np.linspace(10, 40, 12)
    ytheory = 2 / xtheory

    plt.ylim(0,1)
    plt.xlim(0,45)
    plt.show()

    plt.hist(thetas, bins = 18, range=(0,90))
    plt.xlabel('Launch Angle Relative to Surface (degrees)', fontsize = 20)
    plt.ylabel('# Particles', fontsize = 20)
    plt.grid(True, axis = 'y')
    plt.show()

    np.savetxt('Velocity_Distribution.csv', vy, fmt = '%1.3f', delimiter = ',')
    # plot velocity distrubution
    plt.hist(velos, bins=10)
    plt.title('velocity distribution')
    plt.xlabel('vertical velocity (m/s)')
    plt.ylabel('# detections')
    plt.grid(True, axis='y')
    plt.show()

    velocity_slice = []
    for i in range(len(sizes_r)):
       if np.logical_and(30 <= sizes_r[i], sizes_r[i] <= 34):
           velocity_slice.append(velos[i])

    np.savetxt('Velo_dist_30_34_microns.csv', velocity_slice, fmt = '%1.3f', delimiter = ',')


    # plot mass distrubution
    plt.hist(sizes_diam, bins=10)
    plt.title('size distribution')
    plt.xlabel('diameter ($\mu m$)')
    plt.ylabel('# detections')
    plt.grid(True, axis='y')
    plt.show()

# -----------------------------------------------------------------------
# GET THE DATA
# -----------------------------------------------------------------------
skip = [0,1,21,22,23,24]
b = pd.read_csv('master.csv', header=None, skiprows=skip, skip_blank_lines=True)
size_diam = b[8]
bvx = b[9]
bvy = b[10]
bvelos = b[11]
thetas = b[12]
phis = b[13]

plotter(size_diam, bvy, bvx, bvelos, thetas, phis)