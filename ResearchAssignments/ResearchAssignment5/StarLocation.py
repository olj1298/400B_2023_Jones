#olivia jones 5.5.23
#code modified from 400B course by G.Besla

#select galactic center of M31
#selecting particles with similar radius as sun in MW (~8.2kpc)
#selecting particles with similar velocity as sun in MW (~250km/s)
#what snapshots used (Snapshot 800 (∼12 Gyr))
#follow particles over time
#when will stop simulation and relevant to which center

import numpy as np
import astropy
import astropy.units as u
from astropy.constants import G
import astropy.constants as const
from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from OrbitCOM_Soln import orbitCOM
from MassProfile_Soln import MassProfile
import matplotlib.pyplot as plt
import scipy.optimize as so
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
import os
import sys

def RotateFrame(posI,velI): # from in-class lab 7
    """Rotate the position and velocity vectors so that the disk angular momentum is aligned with z axis. 
    Inputs:
        :posI(array of floats): 3D array of positions (x,y,z)
        :velI(array of floats): 3D array of velocities (vx,vy,vz)
    Return:
        ::pos(array of floats): rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        ::vel(array of floats): rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector is in the +z direction"""
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))
    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    # z unit vector
    z_norm = np.array([0, 0, 1])
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2
    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

def solarparticles(path, galaxy, snap, sunparam, suntolr):
    """Indices of solar like particles in galaxies in distance from galactic center and velocity.
        Inputs:
            :path(str): file address data for galaxy
            :galaxy(str): name of galaxy = 'MW' 'M31' 'M33'
            :snap(int): snap number index '0' to '800'
            :parameters(tuple): distance (kpc) and velocity (km/s)
            :sun_tolerances(tuple): tolerances in distance from galactic center and circular speed
        Returns:
            ::indices(numpy array): indices of solar like particles in galaxy"""
    
    sunr, maxzvelocity = sunparam # tuple of sun parameters
    sunrtol, circvtol = suntolr # tuple of sun tolerances
    
    #run mass calculation
    massprof = MassProfile(galaxy, snap)
    calccircspeed = massprof.circularVelocityTotal(sunr).value 
    
    ilbl = '000' + str(snap) # add a string of the filenumber to the value “000”
    ilbl = ilbl[-3:] # remove all but the last 3 digits
    filename = path + galaxy + f'_'+ ilbl + '.txt' # filename of snap being iterated
    
    COMD = CenterOfMass(filename, 2)  #Obtain unitless position and velocity vectors of the COM
    
    COMP = COMD.COM_P(0.1) #define Center of Mass object for COM Position of object at snapshot given
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2]) # define the x, y, and z of velocity vector

    
    # positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 

    # velocities of disk particles relative to COM 
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 

    # Vectors for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T # transposed

    # pos and vel of disk particles post frame of reference rotation so the z-axis aligns with the total angular momentum
    rrot, vrot = RotateFrame(r, v) 
    
    rtoCOM = np.sqrt(rrot[:,0]**2 + rrot[:,1]**2) # x-y distance of particles from COM
    circspeed = np.sqrt(vrot[:,0]**2 + vrot[:,1]**2) # x-y speed of disk particles
    
    # particles with solar like radius in tolerance
    radialindices = np.where(abs((rtoCOM-sunr)/sunr) < sunrtol)
    # particles with solar like circular speed in tolerance
    circvindices = np.where(abs((circspeed-calccircspeed)/calccircspeed) \
                                      < circvtol)
    # particles with z-velocity < 30 km/s in tolerance (van der Marel+ 2012a)
    zindices = np.where(abs(v[:,2]) < maxzvelocity) 
    
    indices = np.intersect1d(radialindices, zindices) # radius and z-velocity cutoffs
    indices = np.intersect1d(indices, circvindices) # circular vel cutoff for MW/M31
    
    return indices

sunr = 8.29 # distance of Sun from MW center (kpc) (van der Marel+ 2012b)
maxzvelocity = 150 # particles with z-velocity < 30 (km/s) (van der Marel+ 2012a)

sunparam = (sunr, maxzvelocity) # storing radial distance of sun from galactic center (kpc) and max velocity (km/s) out of plane
suntolr = (0.2, 0.2) # tolerance values  within 20% percent

snapstart = 0 # start at snap number
snapend = 800 # stop at snap number
increment = 20 # step
#list file paths because for VSCode it needs to know exactly where the file is
MWpath = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'
M31path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'
M33path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'

#run orbit COM for each galaxy
print(f"running orbitCOM for MW")
orbitCOM('MW', snapstart, snapend+1, increment, MWpath)
print(f"running orbitCOM for M31") 
orbitCOM('M31', snapstart, snapend+1, increment, M31path) 
print(f"running orbitCOM for M33")
orbitCOM('M33', snapstart, snapend+1, increment, M33path) 

MWorbit = np.genfromtxt(os.path.abspath('Orbit_MW.txt'), names=True) # create MW COM file
xMW = MWorbit['x'] # x position of MW COM
yMW = MWorbit['y'] # y position of MW COM
zMW = MWorbit['z'] # z position of MW COM

M31orbit = np.genfromtxt(os.path.abspath('Orbit_M31.txt'), names=True) # create M31 COM file
xM31 = M31orbit['x'] # x position of M31 COM
yM31 = M31orbit['y'] # y position of M31 COM
zM31 = M31orbit['z'] # z position of M31 COM

M33orbit = np.genfromtxt(os.path.abspath('Orbit_M33.txt'), names=True) # create M33 COM file
xM33 = M33orbit['x'] # x position of M33 COM
yM33 = M33orbit['y'] # y position of M33 COM
zM33 = M33orbit['z'] # z position of M33 COM


print(f"running MW indices")
MWsunindices = solarparticles(MWpath, 'MW', snapstart, sunparam, suntolr) # solar like particle indices in MW
print(f"running M31 indices")
M31sunindices = solarparticles(M31path, 'M31', snapstart, sunparam, suntolr) # solar like particle indices in M31
print(f"running M33 indices")
M33sunindices = solarparticles(M33path, 'M33', snapstart, sunparam, suntolr) # solar like particle indices in M33

fig, ax = plt.subplots(figsize=(10,10))

for index, i in enumerate(range(snapstart, snapend+1, increment)):
    
    ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
    ilbl = ilbl[-3:] # remove all but the last 3 digits
    
    #hist of radial position of M31 solar like particles
    filename = MWpath + 'MW_' + ilbl + '.txt' # M31 filename including snap
    time, total, data = Read(filename) # read in the datafile
    diskindices = np.where(data['type'] == 2) # get all disk particle indices
    data = data[diskindices] # data for all disk particles
    x = data['x'] - MWorbit['x'][index] # x pos of M31 disk particles relative to COM
    y = data['y'] - MWorbit['y'][index] # y pos of M31 disk particles relative to COM
    z = data['z'] - MWorbit['z'][index] # z pos of M31 disk particles relative to COM
    r = np.sqrt(x**2 + y**2 + z**2) # magnitude of distance vector
    ax.hist(r[MWsunindices], bins=np.linspace(0, 60, 100)) # plot M31 solar like particles
    ax.set_title(f'Number of MW Sunlike Particles at Distance from COM at Snap {i}')
    ax.set_xlabel('Distance from COM [kpc]')
    ax.set_ylabel('Counts')
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 250)    
    
    # save plot as image
    plt.savefig(f'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/MW_migration_histogram_{index}.png') 
    ax.clear() # clear so plots dont stack over snap numbers
    
fig, ax = plt.subplots(figsize=(10,10))

for index, i in enumerate(range(snapstart, snapend+1, increment)):
    
    ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
    ilbl = ilbl[-3:] # remove all but the last 3 digits
    
    #hist of radial position of M31 solar like particles
    filename = M31path + 'M31_' + ilbl + '.txt' # M31 filename including snap
    time, total, data = Read(filename) # read in the datafile
    diskindices = np.where(data['type'] == 2) # get all disk particle indices
    data = data[diskindices] # data for all disk particles
    x = data['x'] - M31orbit['x'][index] # x pos of M31 disk particles relative to COM
    y = data['y'] - M31orbit['y'][index] # y pos of M31 disk particles relative to COM
    z = data['z'] - M31orbit['z'][index] # z pos of M31 disk particles relative to COM
    r = np.sqrt(x**2 + y**2 + z**2) # magnitude of distance vector
    ax.hist(r[M31sunindices], bins=np.linspace(0, 60, 100)) # plot M31 solar like particles
    ax.set_title(f'Number of M31 Sunlike Particles at Distance from COM at Snap {i}')
    ax.set_xlabel('Distance from COM [kpc]')
    ax.set_ylabel('Counts')
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 250)
    
    # save plot as image
    plt.savefig(f'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/M31_migration_histogram_{index}.png') 
    ax.clear() # clear so plots dont stack over snap numbers
    
fig, ax = plt.subplots(figsize=(10,10))

for index, i in enumerate(range(snapstart, snapend+1, increment)):
    ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
    ilbl = ilbl[-3:] # remove all but the last 3 digits
    
    #hist of radial position of M31 solar like particles
    filename = M33path + 'M33_' + ilbl + '.txt' # M31 filename including snap
    time, total, data = Read(filename) # read in the datafile
    diskindices = np.where(data['type'] == 2) # get all disk particle indices
    data = data[diskindices] # data for all disk particles
    x = data['x'] - M33orbit['x'][index] # x pos of M31 disk particles relative to COM
    y = data['y'] - M33orbit['y'][index] # y pos of M31 disk particles relative to COM
    z = data['z'] - M33orbit['z'][index] # z pos of M31 disk particles relative to COM
    r = np.sqrt(x**2 + y**2 + z**2) # magnitude of distance vector
    ax.hist(r[M33sunindices], bins=np.linspace(0, 60, 100)) # plot M31 solar like particles
    ax.set_title(f'Number of M33 Sunlike Particles at Distance from COM at Snap {i}')
    ax.set_xlabel('Distance from COM [kpc]')
    ax.set_ylabel('Counts')
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 250)    
    
    # save plot as image
    plt.savefig(f'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/M33_migration_histogram_{index}.png')
    ax.clear() # clear so plots dont stack over snap numbers
    
# class SolarParticles: #Class to define COM position and velocity properties of a given galaxy #and simulation snapshot

#     def __init__(self,snap,ptype):
#         """Class to calculate the 6-D phase-space position of a galaxy's center of mass using
#         a specified particle type. 
#             Inputs:
#                 :snap(int): snapshot number of file for particle calculation
#                 :ptype(integer): 1, 2, or 3. particle type to use for COM calculations"""
                
#         ilbl = '000' + str(snap) #add a string of the filenumber to the value “000”
#         ilbl = ilbl[-3:] #remove all but the last 3 digits
#         self.filename = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/' + '%s_'%('M31')+ilbl+'.txt'
#         self.time, self.total, self.data = Read(self.filename) #Read in the file
#         self.index = np.where(self.data['type']==ptype) #Create an index to sort by specified ptype
#         #store positions, velocities of only the particles of the given type
#         #Create arrays of all position, velocity components
#         self.x = self.data['x'][self.index]
#         self.y = self.data['y'][self.index]
#         self.z = self.data['z'][self.index]
#         self.vx = self.data['vx'][self.index]
#         self.vy = self.data['vy'][self.index]
#         self.vz = self.data['vz'][self.index]
             
#     def solarselector(self,ptype=2):
#         """This function determines which indices in a snapshot file represent Sun-like stars.
#             Inputs:
#                 :ptype(integer): Particle type, Halo=1, Bulge=2, Disk=3
#             Returns:
#                 ::list of index for particles that have solar properties in file"""
        
#         volDec = 2
#         COM = CenterOfMass(self.filename, ptype) #Obtain unitless position and velocity vectors of the COM
#         pvec = COM.COM_P(0.1,volDec) #define Center of Mass object for COM Position of object at snapshot given
#         velvec = COM.COM_V(pvec[0],pvec[1], pvec[2]).value #define the x, y, and z of velocity vector
#         valpvec = pvec.value #adding units to position vector
#         for i in range(np.size(self.x)): #go through all particles of type specified
#             #Using radius from galactic center defined in (Reid 2014 ApJ 783), (Abuter+2019 A&A 625), Sparke & Gallagher
#             #using velocity (from Marel+2012)
#             sunindex = np.where(7.9 < radius < 9.178 and 239 < velocity < 250)
#             #position and vel vectors
#             calcp = np.array([self.x[sunindex], self.y[sunindex], self.z[sunindex]])
#             calcv = np.array([self.vx[sunindex], self.vy[sunindex], self.vz[sunindex]])
#             #difference of COM to the particles radius and velocity
#             radius = np.linalg.norm(calcp - valpvec) 
#             velocity = np.linalg.norm(calcv - velvec)
#             #lab 7 rotating M31 for xyz plane
#             magdist = np.sqrt(self.pvec[sunindex,0]**2 + self.pvec[sunindex,1]**2 + self.pvec[sunindex,2]**2)
#             magvel = np.sqrt(self.vevec[sunindex,0]**2 + self.velvec[sunindex,1]**2 + self.velvec[sunindex,2]**2)
#             #Rotate vectors so that angular momentum is in z direction
#             rotp, rotv = self.RotateFrame(self.pvec, self.velvec)
#             # Obtain the radial distance in rotated frame
#             rotdist = np.sqrt(rotp[:,0]**2 + rotp[:,1]**2)
#             rotvel = np.sqrt(rotv[:,0]**2 + rotv[:,1]**2 + rotv[:,2]**2)
#             #print(sunindex)
#         return magdist,magvel,rotdist,rotvel
        
#     def RotateFrame(posI,velI): #from Lab7
#         """Rotate the frame of M31 and M33 so edge on and face on components of their componets can be analyzed with z axis. 
#             Inputs:
#                 :posI(float): 3D array of positions
#                 :velI(float): 3D array of velocities
#             Returns:
#                 ::pos positions array in xy plane
#                 ::vel velocities array in z hat angular momentum vector"""
   
#         L = np.sum(np.cross(posI,velI), axis=0) #angular momentum
#         L_norm = L/np.sqrt(np.sum(L**2)) #normalize the vector

#         # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
#         z_norm = np.array([0, 0, 1]) # z unit vector
        
#         # cross product between L and z
#         vv = np.cross(L_norm, z_norm)
#         s = np.sqrt(np.sum(vv**2))
        
#         # dot product between L and z 
#         c = np.dot(L_norm, z_norm)
        
#         # rotation matrix
#         I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#         v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
#         R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

#         # Rotate coordinate system
#         pos = np.dot(R, posI.T).T
#         vel = np.dot(R, velI.T).T
        
#         return pos, vel
    
#     def sunlikestars(path, galaxy, snap, sunr,maxzvelocity,sunrtol,circularspeedtol): #sunparameters, suntolerances):
#         """Get indices solar type stars in the MW, M31, or M33. These stars are defined
#         as ones which are a specified radial distance away from their galactic center and with circular
#         speeds which match the one calculated by the enclosed mass within that radius, within supplied
#         tolerances. There is also a cutoff for maximum out-of-plane velocity (van der Marel+ 2012a).
        
#             Inputs:
#                 :path(str): file path prefix of location of simulation data for galaxy
#                 :galaxy(str):name of galaxy to find indices for ('MW', 'M31', 'M33')
#                 :snap(int): snap number to find sunlike star indices at
#                 :parameters(tuple): tuple containing sun distance from galactic center (kpc) and maximum out-of-plane velocity (km/s)
#                 :sun_tolerances(tuple): fractional allowed tolerances in the radial distance from the galactic center and circular speed
#             Returns:
#                 ::indices: array of indices of sun-like stars in the given galaxy"""
    
#         sunr = float(8.29) #,maxzvelocity
#         maxzvelocity = 100 #sunparameters # unpack tuple of sun parameters
#         sunrtol = float(.1) #, 
#         circularspeedtol = float(.1) #suntolerances # unpack tolerances
#         volDec = 2
#         massprofile = MassProfile(galaxy, snap)
#         calculatedcircularspeed = massprofile.circularVelocityTotal(sunr).value
        
#         ilbl = '000' + str(snap) # add a string of the filenumber to the value “000”
#         ilbl = ilbl[-3:] # remove all but the last 3 digits
#         # compose filename for each snap number being iterated
#         filename = path + f'{galaxy}/' + '%s_'%(galaxy) + ilbl + '.txt' 
        
#         COMD = CenterOfMass(filename, 2) # COM object of galaxy using disk particles
        
#         COMP = COMD.COM_P(0.1,volDec) # get position of COM at this snap
#         COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2]) # velocity of COM at this snap
        
#         # Determine positions of disk particles relative to COM 
#         xD = COMD.x - COMP[0].value 
#         yD = COMD.y - COMP[1].value 
#         zD = COMD.z - COMP[2].value 

#         # Determine velocities of disk particles relative to COM motion
#         vxD = COMD.vx - COMV[0].value 
#         vyD = COMD.vy - COMV[1].value 
#         vzD = COMD.vz - COMV[2].value 

#         # Vectors for r and v 
#         r = np.array([xD,yD,zD]).T # transposed 
#         v = np.array([vxD,vyD,vzD]).T

#         # position and velocity of disk particles after rotating the frame of reference so that
#         # the z-axis lines up with the total galaxy angular momentum
#         rrotated, vrotated = SolarParticles.RotateFrame(r, v) 
        
#         rtoCOM = np.sqrt(rrotated[:,0]**2 + rrotated[:,1]**2) # x-y distance of particles from COM
#         circularspeed = np.sqrt(vrotated[:,0]**2 + vrotated[:,1]**2) # x-y speed of disk particles
        
#         # all particles with sun-like radius within specified tolerance
#         radialindices = np.where(abs((rtoCOM-sunr)/sunr) < sunrtol)
#         # all particles with sun-like circular speed within specified tolerance
#         circularspeedindices = np.where(abs((circularspeed-calculatedcircularspeed)/calculatedcircularspeed) \
#                                         < circularspeedtol)
#         # all particles with z-velocity < 30 km/s within tolerance, condition from (van der Marel+ 2012a)
#         zindices = np.where(abs(v[:,2]) < maxzvelocity) 
        
#         indices = np.intersect1d(radialindices, zindices) # radius and z-velocity cutoffs
#         indices = np.intersect1d(indices, circularspeedindices) # add circular speed cutoff for MW/M31
        
#         return indices
    
# if __name__ == '__main__' :
#     #distances to Galactic Center from the Sun
#     sunr = 8.29 # Distance of Sun from MW center in kpc (van der Marel+ 2012b)
#     maxzvelocity = 100 # all particles with z-velocity < 30 km/s, condition from (van der Marel+ 2012a)
#     volDec = 2.0
#     # tuple storing radial distance of sun from galactic center (kpc) and maximum allowed out-of-plane velocity (km/s)
#     #sunparameters = (sunr, maxzvelocity) 
#     sunrtol = .1
#     circularspeedtol = .1
#     #suntolerances = (0.1, 0.1) # fractional tolerance values for defining sun-like stars
#     #get the radial distance of particles from the COM position (get COM from CenterOfMass.py), circular speeds at that radii
#     #analytic v_circ from enclosed mass is from HW 6, actual tangential velocity may need new codes), and the z-component of
#     #the position (distance from the galactic plane, new code)
#     snapstart = 0 # starting snap number
#     snapend = 100 # ending snap number
#     increment = 1 # increment between snaps for plotting purposes
#     MWpath = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'
#     M31path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'
#     M33path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/'
#     orbitCOM('MW', snapstart, snapend+1, increment)#, MWpath) # MW COM file
#     orbitCOM('M31', snapstart, snapend+1, increment)#, M31path) # M31 COM file
#     orbitCOM('M33', snapstart, snapend+1, increment)#, M33path) # M33 COM file
#     MWorbit = np.genfromtxt(os.path.abspath('Orbit_MW.txt'), names=True) # Read in MW COM file
#     xMW = MWorbit['x'] # x position of MW COM
#     yMW = MWorbit['y'] # y position of MW COM
#     zMW = MWorbit['z'] # z position of MW COM

#     M31orbit = np.genfromtxt(os.path.abspath('Orbit_M31.txt'), names=True) # Read in M31 COM file
#     xM31 = M31orbit['x'] # x position of M31 COM
#     yM31 = M31orbit['y'] # y position of M31 COM
#     zM31 = M31orbit['z'] # z position of M31 COM

#     M33orbit = np.genfromtxt(os.path.abspath('Orbit_M33.txt'), names=True) # Read in M33 COM file
#     xM33 = M33orbit['x'] # x position of M33 COM
#     yM33 = M33orbit['y'] # y position of M33 COM
#     zM33 = M33orbit['z'] # z position of M33 COM

#     MWsunlikeindices = SolarParticles.sunlikestars(MWpath, 'MW',0, str(sunr),str(maxzvelocity),str(sunrtol),str(circularspeedtol)) # parameters,suntolerances) # sun-like particle indices in MW
#     M31sunlikeindices = SolarParticles.sunlikestars(M31path,'M31', 0, str(sunr),str(maxzvelocity),str(sunrtol),str(circularspeedtol)) # sun-like particle indices in M31
#     M33sunlikeindices = SolarParticles.sunlikestars(M33path, 'M33',0, str(sunr),str(maxzvelocity),str(sunrtol),str(circularspeedtol)) # sun-like particle indices in M33
#     # Get positions of all particles at each snapshot relative to the M31 COM. Create frame-by-frame plots at eachsnap number, 
#     # plotting the sun-like particles in a different color.
#     ax = plt.figure(figsize=(10,10)).add_subplot(projection="3d",computed_zorder=False)
#     for index, i in enumerate(range(snapstart, snapend+1, increment)):
#         index = (index + snapstart)*increment # correct index for starting at nonzero snap_start
#         ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
#         ilbl = ilbl[-3:] # remove all but the last 3 digits
        
#         ### Milky Way Plotting
#         filename1 = 'MW' + ilbl + '.txt' # MW filename including snap
#         time, total, data1 = Read(filename1) # read in the datafile
#         diskindices1 = np.where(data1['type'] == 2) # get all disk particle indices
#         data1 = data1[diskindices1] # data for all disk particles
#         x1 = data1['x'] - M31orbit['x'][index] # x position for MW disk particles relative to M31 COM
#         y1 = data1['y'] - M31orbit['y'][index] # y position for MW disk particles relative to M31 COM
#         z1 = data1['z'] - M31orbit['z'][index] # z position for MW disk particles relative to M31 COM
#         ax.scatter(x1, y1, z1, s=0.1, marker='.', c='green', zorder=3) # plot MW disk particles
#         # Optionally, one can include MW sun-like particles in plot. Here, I leave it out.
#         # ax.scatter(x1[MW_sunlike_indices], y1[MW_sunlike_indices],\
#         # z1[MW_sunlike_indices], s=0.1, marker='.', c='black', zorder=5)
        
#         ### M31 Plotting
#         filename2 = 'M31' + ilbl + '.txt' # M31 filename including snap
#         time, total, data2 = Read(filename2) # read in the datafile
#         diskindices2 = np.where(data2['type'] == 2) # get all disk particle indices
#         data2 = data2[diskindices2] # data for all disk particles
#         x2 = data2['x'] - M31orbit['x'][index] # x position for M31 disk particles relative to M31 COM
#         y2 = data2['y'] - M31orbit['y'][index] # y position for M31 disk particles relative to M31 COM
#         z2 = data2['z'] - M31orbit['z'][index] # z position for M31 disk particles relative to M31 COM
#         ax.scatter(x2, y2, z2, s=0.01, marker='.', c='orange') # plot M31 disk particles
#         ax.scatter(x2[M31sunlikeindices], y2[M31sunlikeindices], \
#                 z2[M31sunlikeindices], s=1, c='black') # plot M31 sun-like particles
        
#         ### M33 Plotting
#         filename3 = 'M33' + ilbl + '.txt' # M33 filename including snap
#         time, total, data3 = Read(filename3) # read in the datafile
#         diskindices3 = np.where(data3['type'] == 2)
#         data3 = data3[diskindices3] # data for all disk particles
#         x3 = data3['x'] - M31orbit['x'][index] # x position for M33 disk particles relative to M31 COM
#         y3 = data3['y'] - M31orbit['y'][index] # y position for M33 disk particles relative to M31 COM
#         z3 = data3['z'] - M31orbit['z'][index] # z position for M33 disk particles relative to M31 COM
#         ax.scatter(x3, y3, z3, s=0.01, marker='.', c='red') # plot M33 disk particles
#         ax.scatter(x3[M33sunlikeindices], y3[M33sunlikeindices],\
#                 z3[M33sunlikeindices], s=0.1, c='black') # plot M33 sun-like particles
        
#         ax.grid(True)
#         ax.set_xlim(-50, 50)
#         ax.set_ylim(-50, 50)
#         ax.set_zlim(-50, 50)
#         ax.view_init(30, 30)
        
#         plt.savefig(f'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M31_COM_{index}.png') # save plot as image
#         ax.clear() # clear plot to get the next frame ready
#     fig, ax = plt.subplots(figsize=(10,10))

#     for index, i in enumerate(range(snapstart, snapend+1, increment)):
#         index = (index + snapstart)*increment # correct index for starting at nonzero snap_start
#         ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
#         ilbl = ilbl[-3:] # remove all but the last 3 digits
        
#         ### Plot histograms of radial position of M31 sunlike particles
#         filename = 'M31' + ilbl + '.txt' # M31 filename including snap
#         time, total, data = Read(filename) # read in the datafile
#         diskindicies = np.where(data['type'] == 2) # get all disk particle indices
#         data = data[diskindicies] # data for all disk particles
#         x = data['x'] - M31orbit['x'][index] # x position for M31 disk particles relative to M31 COM
#         y = data['y'] - M31orbit['y'][index] # y position for M31 disk particles relative to M31 COM
#         z = data['z'] - M31orbit['z'][index] # z position for M31 disk particles relative to M31 COM
#         r = np.sqrt(x**2 + y**2 + z**2)
#         print(np.max(r))
#         ax.hist(r[M31sunlikeindices], bins=np.linspace(0, 60, 100)) # plot M31 sun-like particles
#         ax.set_title(f'Number of M31 Sunlike Particles at Given Distances from M31 COM (kpc) at Snap {i}')
#         ax.set_xlabel('Distance from COM [kpc]')
#         ax.set_ylabel('Counts')
#         ax.set_xlim(0, 60)
#         ax.set_ylim(0, 200)    
#         plt.savefig(f'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M31_migration_histogram_{index}.png') # save plot as image
#         ax.clear() # clear plot to get the next frame ready