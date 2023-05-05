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
import astropy.constants as const
from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from OrbitCOM_Soln import orbitCOM
from MassProfile_Soln import MassProfile
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import os

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

#user changeable variables depending on simulation
sunr = 8.29 # distance of Sun from MW center (kpc) (van der Marel+ 2012b)
maxzvelocity = 30 # particles with z hat velocity < 30 (km/s) (van der Marel+ 2012a)
sunparam = (sunr, maxzvelocity) # storing radial distance of sun from galactic center (kpc) and max velocity (km/s) out of plane
suntolr = (0.2, 0.2) # tolerance values  within 20% percent

snapstart = 0 # start at snap number
snapend = 800 # stop at snap number
increment = 20 # step
#list file paths because for VSCode it needs to know exactly where the file is
MWpath = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment7/'
M31path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment7/'
M33path = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment7/'

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


print(f"running galaxies indices, check folder for graphs to update")
MWsunindices = solarparticles(MWpath, 'MW', snapstart, sunparam, suntolr) # solar like particle indices in MW
M31sunindices = solarparticles(M31path, 'M31', snapstart, sunparam, suntolr) # solar like particle indices in M31
M33sunindices = solarparticles(M33path, 'M33', snapstart, sunparam, suntolr) # solar like particle indices in M33

fig, ax = plt.subplots(figsize=(10,10))

for index, i in enumerate(range(snapstart, snapend+1, increment)):
    
    ilbl = '000' + str(i) # add a string of the filenumber to the value “000”
    ilbl = ilbl[-3:] # remove all but the last 3 digits
    
    #radial position of M31 solar like particles
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