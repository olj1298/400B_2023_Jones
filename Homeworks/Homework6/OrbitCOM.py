#Homework 6 Template
#G. Besla & R. Li
#o.jones 2.23.23
import numpy as np #import modules
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt # import plotting modules
import matplotlib
import sys
import os
from ReadFile import Read #my modules
from CenterOfMass import CenterOfMass
#Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying by how much to decrease RMAX instead of a factor of 2

def OrbitCOM(galaxy,start,end,n):
    """Calculate COM position and velocity as a function of time over defined snapshots.
        Inputs:
            :galaxy(string): galaxy name 'MW'
            :start(float): number of first snapshot to read
            :end(float): number of last snapshot to be read
            :n(integer): intervals over which we return to COM
        Returns: 
            ::save orbits output as file"""

    fileout = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/'+'Orbit_'+galaxy+'_'+str(start)+'_'+str(end)+'.txt' #savefile for orbit calculated
    #print(f"fileout: {fileout}")
    delta = 0.1 #tolerance
    volDec = 2 #amount rmax is decreased in CenterOfMass function
    snapids = np.arange(start,end+1,n) #create array
    arraysize = np.size(snapids) #for array size check
    orbit = np.zeros([len(snapids),7]) #array stores x,y,z,vx,vy,vz of COM of each object
    
    if arraysize == 0: #check if array is empty and code contains error. Break if so
        print(f"array is zero, check code.")
        sys.exit()
    
    for  i,snapid in enumerate(snapids): #loop for all snapshots
        
        ilbl = '000'+str(snapid)
        ilbl = ilbl[-3:] #keep last 3 digits
        filename = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/'+(galaxy)+'_'+str(ilbl)+'.txt' #save to directory address
        #print(f"filename is: {filename}")
        COM = CenterOfMass(filename,2) #Instance of CenterOfMass with 2 = disk particles

        if galaxy == 'M33': #condition for M33
            delta = 0.1
            volDec = 4 #M33 tidally stripped towards end of simulation
            COM_p = COM.COM_P(delta,volDec)
            COM_v = COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
        else:
            delta = 0.1
            volDec = 2
            COM_p = COM.COM_P(delta,volDec)
            COM_v = COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
            
        time = COM.time/1000 #time, pos, vel in ith element of the orbit array, without units 
        orbit[i][0] = time.value
        orbit[i][1] = COM_p[0].value
        orbit[i][2] = COM_p[1].value
        orbit[i][3] = COM_p[2].value
        orbit[i][4] = COM_v[0].value
        orbit[i][5] = COM_v[1].value
        orbit[i][6] = COM_v[2].value
        
        print(snapid) #test
        
    #write the data to a file
    #we do this because we don't want to have to repeat this process 
    #this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    #return

#Recover the orbits and generate the COM files for each galaxy read in 800 snapshots in intervals of n=5
#Note: This might take a little while - test your code with a smaller number of snapshots first!
OrbitCOM('MW',0,800,5)
OrbitCOM('M31',0,800,5)
OrbitCOM('M33',0,800,5)
#Read in the data files for the orbits of each galaxy that you just created
#headers: t, x, y, z, vx, vy, vz
#using np.genfromtxt

dataMW = np.genfromtxt(os.path.abspath("Orbit_MW_0_800.txt"),dtype=None,names=True,skip_header=0)
dataM31 = np.genfromtxt(os.path.abspath("Orbit_M31_0_800.txt"),dtype=None,names=True,skip_header=0)
dataM33 = np.genfromtxt(os.path.abspath("Orbit_M33_0_800.txt"),dtype=None,names=True,skip_header=0)

#function to compute the magnitude of the difference between two vectors 
#You can use this function to return both the relative position and relative velocity for two 
#galaxies over the entire orbit  

def VectorDifference(vector1,vector2):
    """Computes magnitude difference between 2 vectors.  
    Inputs
        :vector1(narray): coords. of first vector
        :vector2(narray): coords. of second vector
    Returns:
        ::magnitude difference 2 vectors"""
        
    diff = vector1-vector2 #difference of 2 vectors
    magdifference = 0 #begin counter
    for x in diff: #loop of difference vector
        magdifference += x**2 #sum of squares
    magdiff= np.sqrt(magdifference) #square root of sum of squares
    
    return magdiff

#Determine the magnitude of the relative position and velocities 

#MW mag of position and vel
MWpositions= np.zeros([len(dataMW['x']),3])
for i in range(0,len(dataMW['x'])):
    MWpositions[i] = [dataMW['x'][i],dataMW['y'][i],dataMW['z'][i]]

MWvelos= np.zeros([len(dataMW['vx']),3])
for i in range(0,len(dataMW['vx'])):
    MWvelos[i] = [dataMW['vx'][i],dataMW['vy'][i],dataMW['vz'][i]]

#M31 mag of position and vel
M31positions= np.zeros([len(dataM31['x']),3])
for i in range(0,len(dataM31['x'])):
    M31positions[i] = [dataM31['x'][i],dataM31['y'][i],dataM31['z'][i]]

M31velos= np.zeros([len(dataM31['vx']),3])
for i in range(0,len(dataM31['vx'])):
    M31velos[i] = [dataM31['vx'][i],dataM31['vy'][i],dataM31['vz'][i]]

#M33 mag of position and vel
M33positions= np.zeros([len(dataM33['x']),3])
for i in range(0,len(dataM33['x'])):
    M33positions[i] = [dataM33['x'][i],dataM33['y'][i],dataM33['z'][i]]

M33velos= np.zeros([len(dataM33['vx']),3])
for i in range(0,len(dataM33['vx'])):
    M33velos[i] = [dataM33['vx'][i],dataM33['vy'][i],dataM33['vz'][i]]

#Difference of MW and M31
differencepositionsMWM31 = np.array([])
differencevelosMWM31 = np.array([])
for i in range(0,len(dataMW['x'])):
    differencepositionsMWM31 = np.append(differencepositionsMWM31,VectorDifference(MWpositions[i],M31positions[i]))
    differencevelosMWM31 = np.append(differencevelosMWM31,VectorDifference(MWvelos[i], M31velos[i]))

#Difference of M33 and M31
differencepositionsM33M31 = np.array([]) #initialize arrays
differencevelosM33M31 = np.array([])
for i in range(0,len(dataM33['x'])):
    differencepositionsM33M31 = np.append(differencepositionsM33M31,VectorDifference(M33positions[i],M31positions[i]))
    differencevelosM33M31 = np.append(differencevelosM33M31,VectorDifference(M33velos[i],M31velos[i]))

#Plot the Orbit of the galaxies

#M33andM31
plt.plot(dataMW['t'],differencepositionsM33M31,color='red',linewidth=2,label='M33-M31')
plt.xlabel('time (Gyr)')
plt.ylabel('Distance (Kpc)')
plt.title('Orbit of M31-M33')
plt.legend()
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/OrbitofM31andM33')

plt.plot(dataMW['t'],differencepositionsMWM31,color='red',linewidth=2,label='MW-M31')
plt.xlabel('time (Gyr)')
plt.ylabel('Distance (Kpc)')
plt.title('Orbit of MW-M31')
plt.legend()
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/OrbitofMWandM31')
 
#Plot the orbital velocities of the galaxies 

plt.plot(dataMW['t'],differencevelosM33M31,color='red',linewidth=2,label='M33-M31')
plt.xlabel('time (Gyr)')
plt.ylabel('Velocity (km/s)')
plt.title('Relative Velocity of M31-M33')
plt.legend()
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/VeloofM31andM33')

plt.plot(dataMW['t'],differencevelosMWM31,color='red',linewidth=2,label='MW-M31')
plt.xlabel('time (Gyr)')
plt.ylabel('Velocity (km/s)')
plt.title('Relative Velocity of MW-M31')
plt.legend()
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework6/VeloofMWandM31')

"""Questions 1)Two encounters before they merge.
    2)When the galaxies approach each other, the velocity increases. 
    3)M31 and MW around 6 Gyr from now.
    bonus)M33 orbit decays by 10 kpc/Gyr at 6 Gyr. If constant, it would take 7.5 Gyr for M33 to 
    merge with the combined MW and M31."""