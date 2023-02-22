# Homework 6 Template
# G. Besla & R. Li
#import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
import sys
# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass import CenterOfMass

def OrbitCOM(galaxy,start,end,n):
    """Function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
        Inputs:
            :galaxy(string): galaxy name 'MW'
            :start(float): number of first snapshot to read
            :end(float): number of last snapshot to be read
            :n(integer): intervals over which we return to COM
        Returns: 
            ::save orbits output as file"""

    #set tolerance and VolDec for calculating COM_P in CenterOfMass
    #for M33 that is stripped more, use different values for VolDec
    #generate the snapshot id sequence 
    #it is always a good idea to also check if the input is eligible (not required)
    #initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    fileout= 'Orbit_'+galaxy+'_'+str(start)+'_'+str(end)+'.txt'
    
    delta=.1
    volDec=2
           
    snap_ids = np.arange(start,end+1,n) 
    arraysize = np.size(snap_ids)
    
    if arraysize == 0: #check if array is empty and code contains error. Break if so
        print(f"array is zero, check code.")
        sys.exit()

    orbit = np.zeros(len(snap_ids),7) #array stores x,y,z,vx,vy,vz of COM of each object
    
    for  i,snap_id in enumerate(snap_id):# loop over files
        
        ilbl = '000'+str(snap_id)
        ilbl = ilbl[-3:] #remove all but last 3 digits
        filename= 'C/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6'+galaxy+'_' + str(ilbl) +'.txt'
        
        COM = CenterOfMass(filename, 2) # Initialize an instance of CenterOfMass class, using disk particles

        if galaxy == 'M33': #condition for M33
            volDec = 4 #M33 tidally stripped towards end of simulation
        else:
            COM_p = COM.COM_P(delta,volDec)
            COM_v= COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
            
        time = COM.time / 1000 # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        orbit[i][0] =time.value
        orbit[i][1] =COM_p[0].value
        orbit[i][2] =COM_p[1].value
        orbit[i][3] =COM_p[2].value
        orbit[i][4] =COM_v[0].value
        orbit[i][5] =COM_v[1].value
        orbit[i][6] =COM_v[2].value
        
        #a[i] = var1, *tuple(array1)
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    return

# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first!
# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

data_MW = np.genfromtxt("Orbit_MW_0_800.txt",dtype=None,names=True,skip_header=0)
data_M31=np.genfromtxt("Orbit_M31_0_800.txt",dtype=None,names=True,skip_header=0)
data_M33=np.genfromtxt("Orbit_M33_0_800.txt",dtype=None,names=True,skip_header=0)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def VectorDifference(vector1,vector2):
    """Computes difference between 2 vectors and then find magnitude of the difference vector
    Inputs
        :vector1():
        :vector2():
    Returns:
        ::magnitudedifference(float) mag of diff of 2 vectors that were pass in"""
        
    dif = vector1-vector2 # difference of 2 vectors
    magnitude_difference = 0 # initialize the magnitude
    for x in dif: # loop through difference vector
        magnitude_difference+= x**2 # sum the squares
    magnitude_difference= np.sqrt(magnitude_difference) # square root of the sum of squares
    
    return magnitude_difference

# Determine the magnitude of the relative position and velocities 
# of MW and M31
# of M33 and M31

# Determine the magnitude of the relative position and velocities 

MW_positions= np.zeros([len(data_MW['x']),3])
for i in range(0,len(data_MW['x'])):
    MW_positions[i] = [data_MW['x'][i],data_MW['y'][i],data_MW['z'][i]]

M31_positions= np.zeros([len(data_M31['x']),3])
for i in range(0,len(data_M31['x'])):
    M31_positions[i] = [data_M31['x'][i],data_M31['y'][i],data_M31['z'][i]]


M33_positions= np.zeros([len(data_M33['x']),3])
for i in range(0,len(data_M33['x'])):
    M33_positions[i] = [data_M33['x'][i],data_M33['y'][i],data_M33['z'][i]]





MW_velos= np.zeros([len(data_MW['vx']),3])
for i in range(0,len(data_MW['vx'])):
    MW_velos[i] = [data_MW['vx'][i],data_MW['vy'][i],data_MW['vz'][i]]

M31_velos= np.zeros([len(data_M31['vx']),3])
for i in range(0,len(data_M31['vx'])):
    M31_velos[i] = [data_M31['vx'][i],data_M31['vy'][i],data_M31['vz'][i]]


M33_velos= np.zeros([len(data_M33['vx']),3])
for i in range(0,len(data_M33['vx'])):
    M33_velos[i] = [data_M33['vx'][i],data_M33['vy'][i],data_M33['vz'][i]]

# of MW and M31
difference_positions_MW_M31 = np.array([])
difference_velos_MW_M31 = np.array([])
for i in range(0,len(data_MW['x'])):
    difference_positions_MW_M31 = np.append(difference_positions_MW_M31,VectorDifference(MW_positions[i], M31_positions[i]))
    difference_velos_MW_M31 = np.append(difference_velos_MW_M31,VectorDifference(MW_velos[i], M31_velos[i]))




# of M33 and M31
difference_positions_M33_M31 = np.array([]) #initialize arrays
difference_velos_M33_M31 = np.array([])
for i in range(0,len(data_M33['x'])):
    difference_positions_M33_M31 = np.append(difference_positions_M33_M31,VectorDifference(M33_positions[i], M31_positions[i]))
    difference_velos_M33_M31 = np.append(difference_velos_M33_M31,VectorDifference(M33_velos[i], M31_velos[i]))


# Plot the Orbit of the galaxies 
#################################
#fig,ax = plt.subplots(figsize = (10,10))
##ax.semilogy(r,MWH,label = 'Halo Profile')
#ax.semilogy(r,MWD,label = 'Disk Profile')
#ax.semilogy(r,MWB,label = 'Bulge Profile')
#ax.semilogy(r,MWT,label = 'Total Profile')
#ax.semilogy(r,Hern,linestyle = 'dashed',label = 'Hernquist Profile')
#legend = ax.legend() #create legend
#ax.set(title = 'MW Mass Profiles - Scale fit 62 kpc', xlabel = 'Radius (kpc)', ylabel = 'Log(Mass Enclosed ($M_{\odot}$))') #label graph axes
#print(f"The best fit scale length for the MW is: {scaleMW*u.kpc}")
#plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/MWMassProfile.png')

plt.plot(data_MW['t'],difference_positions_M33_M31, color = 'red', linewidth =2, label = 'M33-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Distance (Kpc)')
# plt.title('Orbit of M31-M33')
# plt.legend()
# plt.savefig('Orbit_M31_M33')
# 
# plt.plot(data_MW['t'],difference_positions_MW_M31, color = 'red', linewidth =2, label = 'MW-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Distance (Kpc)')
# plt.title('Orbit of MW-M31')
# plt.legend()
# plt.savefig('Orbit_MW_M31')
# 
# 
# 
# 
# # Plot the orbital velocities of the galaxies 
# #################################
# plt.plot(data_MW['t'],difference_velos_M33_M31, color = 'red', linewidth =2, label = 'M33-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Velocity (km/s)')
# plt.title('Relative Velocity of M31-M33')
# plt.legend()
# plt.savefig('Velo_M31_M33')
# 
# plt.plot(data_MW['t'],difference_velos_MW_M31, color = 'red', linewidth =2, label = 'MW-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Velocity (km/s)')
# plt.title('Relative Velocity of MW-M31')
# plt.legend()
# plt.savefig('Velo_MW_M31')
# =============================================================================

'''Questions 1) There will be 2 close encounters between MW and M33 before they merger in the third encounter
    2) As the galaxies get closer, the velocity increases. Essentially, velocity appears inversely related to
    separation.
3) MW and M31 fully merge just after 6 Gyr from now
bonus) M33 appears to be decaying by ~ 10 kpc / Gyr just after 6 Gyr. If this continued, it would take ~7.5 Gyrs for M33 to 
    then merge with the MW-M31 remnant.

# Plot the orbital velocities of the galaxies 
#################################'''