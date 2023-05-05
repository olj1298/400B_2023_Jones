# Homework 6 Solution
# G. Besla & R. Li & H. Foote
# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
import os
# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass

def orbitCOM(galaxy, start, end, n, path):
    """Loops over all the desired snapshots to compute the COM pos and vel as a function of time.
        Inputs:
            :galaxy(str): name of galaxy 'MW', 'M31', 'M33'
            :start(int): number of first snapshot to be analyzed 
            :end(int): number of last snapshot to be analyzed
            :n(int): increment between each successive snapshot number to be analyzed
            :path(str): relative path between code directory and snap files
        Returns: 
            ::datafile 'Orbit_[galaxy].txt' containing time series with array of position components and velocity components. 
                Saved to directory, not returned"""
    
    # compose the filename for output
    fileout = f"Orbit_{galaxy}.txt"

    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1 # tolerance for iterative COM calculation
    volDec = 2 # factor to shrink sphere for iterative COM calculation
    if galaxy == 'M33': # stricter factor to shrink sphere for M33 due to severe tidal stripping
        volDec = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end+1, n) # array of snap numbers (include endpoint)
    if len(snap_ids) < 1: # check if array is valid
        raise ValueError('Bad array of snapshots: make sure that start < end and n is positive')
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((len(snap_ids), 7)) # empty array to store parameters of orbit
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids):  # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id) # add a string of the filenumber to the value “000”
        ilbl = ilbl[-3:] # remove all but the last 3 digits
        # compose filename for each snap number being iterated
        filename = path + f'%s_'%(galaxy) + ilbl + '.txt' #'{galaxy}/' + 
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, ptype=2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_p = COM.COM_P(delta, volDec) # galaxy position COM
        COM_v = COM.COM_V(COM_p[0], COM_p[1], COM_p[2]) # galaxy COM velocity
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value/1000, *tuple(COM_p.value), *tuple(COM_v.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

# if __name__ == 'main':

#     # Recover the orbits and generate the COM files for each galaxy
#     # read in 800 snapshots in intervals of n=5
#     # Note: This might take a little while - test your code with a smaller number of snapshots first! 
#     orbitCOM('MW', 0, 800, 5) # MW COM file
#     orbitCOM('M31', 0, 800, 5) # M31 COM file
#     orbitCOM('M33', 0, 800, 5) # M33 COM file

#     # Read in the data files for the orbits of each galaxy that you just created
#     # headers:  t, x, y, z, vx, vy, vz
#     # using np.genfromtxt
#     MW_orbit = np.genfromtxt(os.path.abspath('Orbit_MW.txt'), names=True) # Read in MW COM file
#     M31_orbit = np.genfromtxt(os.path.abspath('Orbit_M31.txt'), names=True) # Read in M31 COM file
#     M33_orbit = np.genfromtxt(os.path.abspath('Orbit_M33.txt'), names=True) # Read in M33 COM file

#     # function to compute the magnitude of the difference between two vectors 
#     # You can use this function to return both the relative position and relative velocity for two 
#     # galaxies over the entire orbit  
#     def vector_magnitude(v1, v2):
#         """Determines the magnitude of the difference in two vectors. If the vectors are 2D
#         arrays of the same length, then the function will compute the magnitude of the difference vector
#         at each row in the input.
#             Inputs:
#                 :v1(np.array): first vector to be analyzed
#                 :v2(np.array): second vector to be analyzed
#             Returns: 
#                 ::'float' or 'numpy array' if the vectors are 1D, float of the magnitude of their difference.
#                     if the vectors are 2D, array of difference vector magnitudes."""
                    
#         if v1.shape != v2.shape: # check if vectors have same dimension and length
#             raise ValueError("Input vectors must have same shape.")
#         return np.linalg.norm((v1-v2), axis=0) # return magnitude of each row (single row if 1D vectors)

#     # Determine the magnitude of the relative position and velocities 
#     MW_position = np.array((MW_orbit['x'], MW_orbit['y'], MW_orbit['z'])) # Position vector of MW COM
#     M31_position = np.array((M31_orbit['x'], M31_orbit['y'], M31_orbit['z'])) # Position vector of M31 COM
#     M33_position = np.array((M33_orbit['x'], M33_orbit['y'], M33_orbit['z'])) # Position vector of M33 COM
#     MW_velocity = np.array((MW_orbit['vx'], MW_orbit['vy'], MW_orbit['vz'])) # Velocity vector of MW COM
#     M31_velocity = np.array((M31_orbit['vx'], M31_orbit['vy'], M31_orbit['vz'])) # Velocity vector of M31 COM
#     M33_velocity = np.array((M33_orbit['vx'], M33_orbit['vy'], M33_orbit['vz'])) # Velocity vector of M33 COM

#     # of MW and M31
#     MW_M31_separation_distance = vector_magnitude(MW_position, M31_position) # Relative separation between MW and M31
#     MW_M31_relative_velocity_magnitude = vector_magnitude(MW_velocity, M31_velocity) # Relative speed between MW and M31

#     # of M33 and M31
#     M33_M31_separation_distance = vector_magnitude(M33_position, M31_position) # Relative separation between M33 and M31
#     M33_M31_relative_velocity_magnitude = vector_magnitude(M33_velocity, M31_velocity) # Relative speed between M33 and M31

#     # Plot the Orbit of the galaxies 
#     #################################
#     fig, ax = plt.subplots(1, 2, figsize=(15, 5)) # 2 plots, one for each galaxy pair
#     ax[0].plot(MW_orbit['t'], MW_M31_separation_distance, lw=2) # MW-M31 separation distance
#     ax[0].set_title('MW-M31 Separation Distance (kpc) vs Time (Gyr)')
#     ax[0].set_xlabel('Time (Gyr)')
#     ax[0].set_ylabel('Separation Distance (kpc)')
#     ax[1].plot(MW_orbit['t'], M33_M31_separation_distance, c='red', lw=2) # M31-M33 separation distance
#     ax[1].set_title('M31-M33 Separation Distance (kpc) vs Time (Gyr)')
#     ax[1].set_xlabel('Time (Gyr)')
#     ax[1].set_ylabel('Separation Distance (kpc)')

#     # Plot the orbital velocities of the galaxies 
#     #################################
#     fig, ax = plt.subplots(1, 2, figsize=(15, 5)) # 2 plots, one for each galaxy pair
#     ax[0].plot(MW_orbit['t'], MW_M31_relative_velocity_magnitude, lw=2) # MW-M31 relative speed
#     ax[0].set_title('MW-M31 Relative Velocity Magnitude (km/s) vs Time (Gyr)')
#     ax[0].set_xlabel('Time (Gyr)')
#     ax[0].set_ylabel('Relative Velocity Magnitude (km/s)')
#     ax[1].plot(MW_orbit['t'], M33_M31_relative_velocity_magnitude, c='red', lw=2) # M31-M33 relative speed
#     ax[1].set_title('M31-M33 Relative Velocity Magnitude (km/s) vs Time (Gyr)')
#     ax[1].set_xlabel('Time (Gyr)')
#     ax[1].set_ylabel('Relative Velocity Magnitude (km/s)')