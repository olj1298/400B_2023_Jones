# In Class Lab 7 Template
# G. Besla
# with code from R. Hoffman and E. Patel

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import os
import sys
# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# # This Lab:
# A) Using Contour plots to identify substructure within the stellar disk of M31.
# B) Rotating the disk so that we are seeing it edge on
# C) Create plots to examine the kinematics of the disk

import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """Create a density contour plot.
        Inputs:
            :xdata(array):
            :ydata(array):
            :nbins_x(int): Number of bins along x dimension
            :nbins_y(int): Number of bins along y dimension
            :ax(plt.Axes): If supplied, plot the contour to this axis. Otherwise, open a new figure
            :contour_kwargs(dict): kwargs to be passed to pyplot.contour()
        Returns:
            ::density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)e.g.:
            density_contour(xD, yD, 80, 80, ax=ax,colors=['red','orange', 'yellow', 'orange', 'yellow'])"""

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))

    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.68','0.95', '0.99'][::-1]
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour


# Use the CenterOfMass code to compute the positions and velocities of all particles in M31's disk relative to its center of mass position and motion.

# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass(os.path.abspath("M31_000.txt"),2)

# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 

# total velocity 
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors for r and v 
r = np.array([xD,yD,zD]).T # transposed 
v = np.array([vxD,vyD,vzD]).T

# # Part A:
# Create plot of M31's disk density, using 2D Histograms 
# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 using a 2D historgram
plt.hist2d(xD,yD, bins=200, norm=LogNorm(), cmap='plasma')

plt.colorbar()

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added

density_contour(xD,yD,80,80,ax=ax,colors=['magenta','cyan','yellow','orange']) #change colors based on number of sigma values

# Add axis labels
plt.xlabel(' ', fontsize=22)
plt.ylabel(' ', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save to a file
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/Lab7_M31Disk.png')

#Part B
# Utilize the below code to rotate the M31 disk and plot it edge on and face on. What is the sense of rotation of M31 ? 

def RotateFrame(posI,velI):
    """Function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
        Inputs:
            :posI(float): 3D array of positions (x,y,z)
            :velI(float): 3D array of velocities (vx,vy,vz)
        Returns:
            ::pos: rotated 3D array of positions (x,y,z) such that disk is in the XY plane
            ::vel: rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector is in the +z direction"""
   
    L = np.sum(np.cross(posI,velI), axis=0) # compute the angular momentum
    L_norm = L/np.sqrt(np.sum(L**2)) # normalize the vector

    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    z_norm = np.array([0, 0, 1]) # z unit vector
    
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

# determine the rotated velocity vectors
rn,vn = RotateFrame(r,v)

# Rotated M31 Disk - EDGE ON
# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 , 2D histogram
plt.hist2d(rn[:,0], rn[:,2], bins=200, norm=LogNorm(), cmap='magma')
plt.colorbar()

# Add axis labels
plt.xlabel(' ', fontsize=22)
plt.ylabel(' ', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save to a file
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/Lab7_EdgeOn_Density.png')

# Rotated M31 Disk - FACE ON
# M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for M31 
plt.hist2d(rn[:,0], rn[:,1], bins=200, norm=LogNorm(), cmap='magma')
plt.colorbar()

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
density_contour(rn[:,0], rn[:,1],80,80,ax=ax,colors=['magenta','cyan','yellow','orange']) #change colors based on number of sigma values
# Add axis labels
plt.xlabel('  ', fontsize=22)
plt.ylabel('  ', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save to a file 
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/Lab7_FaceOn_Density.png')

# Part C
# a) Create a scatter plot of the edge on disk particles, weighted by velocity.
# Plot velocity weighted EDGE ON DISK

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# plot position of disk particles color coded by velocity along the 3rd axis
# plt.scatter(pos1, pos2, c=vel1)
plt.scatter(rn[:,0], rn[:,2], c=vn[:,1]) 

#colorbar
cbar = plt.colorbar()
cbar.set_label('  ', size=22)
density_contour(rn[:,0], rn[:,2],80,80,ax=ax,colors=['magenta','cyan','yellow','orange']) #change colors based on number of sigma values
# Add axis labels
plt.xlabel('  ', fontsize=22)
plt.ylabel('  ', fontsize=22)

# calculate the 2D density of the data given
#counts,xbins,ybins=np.histogram2d(xD[index],yD[index],bins=100,normed=LogNorm())

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

# Save file
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/Lab7_EdgeOn_Vel.png')

# b) Create a phase diagram for the edge on disk (Position vs. Velocity) using a 2D Histogram.
# c) Use the MassProfile Code to overplot the expected circular velocity from the mass profile. 
M31 = MassProfile("M31",0)
R = np.arange(0.01,40,0.1)
Vcirc = M31.CircularVelocityTotal(R)

# Make a phase diagram
# MW Disk Velocity Field edge on.

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot 2D Histogram one component of  Pos vs Vel 
plt.hist2d(rn[:,0], vn[:,1], bins=200, norm=LogNorm(), cmap='magma')

# Overplot Circular Velocity from the MassProfile Code
density_contour(rn[:,0], vn[:,1],80,80,ax=ax,colors=['magenta','cyan','yellow','orange']) #change colors based on number of sigma values
plt.plot(R,Vcirc, color='red')
plt.plot(-R,-Vcirc,color='red')

# Add axis labels
plt.xlabel(' ', fontsize=22)
plt.ylabel(' ', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save file
plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/Lab7_PhaseDiagram.png')