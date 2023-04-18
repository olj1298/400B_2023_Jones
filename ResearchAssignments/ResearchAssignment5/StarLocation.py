#olivia jones 4.17.23
#code modified from 400B course by G.Besla
"""Fate of Stars at Sun's Location in the M31 Disk.
Using the known values for our Sun's radius and velocity in the Milky Way provided in Martinez-Barbosa+ 2015, 
we can extend this to similar particles that exist in M31. I hope to come to a conclusion on how to select 
sun like particles. I plan to plot circular velocity vs. radius from galactic center."""

#select galactic center of M31
#selecting particles with similar radius as sun in MW (~8.2kpc)
#selecting particles with similar velocity as sun in MW (~250km/s)
#what snapshots used (Snapshot 800 (∼12 Gyr))
#follow particles over time
#when will stop simulation and relevant to which center

import numpy as np
import astropy.units as u
import astropy.constants as const
from ReadFile import Read
from pathlib import Path
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
import scipy.optimize as so
# import plotting modules
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# my modules
from MassProfile import *
from Lab7 import *
import os
import sys

#1. Water Maser Distance for the Sun :  R_o = 8.34 kpc (Reid 2014 ApJ 783) 
#2. GRAVITY Collaboration Distance for the Sun:  R_o = 8.178 kpc (Abuter+2019 A&A 625)
#3. Value for Distance to Sun listed in Sparke & Gallagher : R_o = 7.9 kpc 
#distances to Galactic Center from the Sun
RoWM = 8.34 * (u.kpc) #kpc
RoG = 9.178 * (u.kpc) #kpc
RoSG = 7.9 * (u.kpc) #kpc

class SolarParticles: #Class to define COM position and velocity properties of a given galaxy #and simulation snapshot

    def __init__(self,snap,ptype):
        """Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            Inputs:
                :snap(int): snapshot number of file for particle calculation
                :ptype(integer): 1, 2, or 3. particle type to use for COM calculations"""
        ilbl = '000' + str(snap) #add a string of the filenumber to the value “000”
        ilbl = ilbl[-3:] #remove all but the last 3 digits
        self.filename = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/' + '%s_'%('M31')+ilbl+'.txt'
        self.time, self.total, self.data = Read(self.filename) #Read in the file
        self.index = np.where(self.data['type']==ptype) #Create an index to sort by specified ptype
        # store the mass, positions, velocities of only the particles of the given type
        #self.m = self.data['m'][self.index] #data in m column from file.
        #Create arrays of all position, velocity components
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.pvec = np.array([self.x,self.y,self.z]).T #transpose array
        self.velvec = np.array([self.vx,self.vy,self.vz]).T #transpose array
        
    def solarselector(self,ptype=2):
        """This function determines which indices in a snapshot file represent Sun-like stars.
            Inputs:
                :ptype(integer): Particle type, Halo=1, Bulge=2, Disk=3
            Returns:
                ::list of index for particles that have solar properties in file"""
        
        COM = CenterOfMass(self.filename, ptype) #Obtain unitless position and velocity vectors of the COM
        pvec = COM.COM_P(0.1) #define Center of Mass object for COM Position of object at snapshot given
        velvec = COM.COM_V(pvec[0],pvec[1], pvec[2]).value #define the x, y, and z of velocity vector
        valpvec = pvec.value #adding units to position vector
        for i in range(np.size(self.x)): #go through all particles of type specified
            #Using radius from galactic center defined in (Reid 2014 ApJ 783), (Abuter+2019 A&A 625), Sparke & Gallagher
            #using velocity (from Marel+2012)
            sunindex = np.where(7.9 < radius < 9.178 and 239 < velocity < 250)
            #position and vel vectors
            calcp = np.array([self.x[sunindex], self.y[sunindex], self.z[sunindex]])
            calcv = np.array([self.vx[sunindex], self.vy[sunindex], self.vz[sunindex]])
            #difference of COM to the particles radius and velocity
            radius = np.linalg.norm(calcp - valpvec) 
            velocity = np.linalg.norm(calcv - velvec)
            #lab 7 rotating M31 for xyz plane
            magdist = np.sqrt(self.pvec[sunindex,0]**2 + self.pvec[sunindex,1]**2 + self.pvec[sunindex,2]**2)
            magvel = np.sqrt(self.vevec[sunindex,0]**2 + self.velvec[sunindex,1]**2 + self.velvec[sunindex,2]**2)
            #Rotate vectors so that angular momentum is in z direction
            rotp, rotv = self.RotateFrame(self.pvec, self.velvec)
            # Obtain the radial distance in rotated frame
            rotdist = np.sqrt(rotp[:,0]**2 + rotp[:,1]**2)
            rotvel = np.sqrt(rotv[:,0]**2 + rotv[:,1]**2 + rotv[:,2]**2)
            print(sunindex)
          
if __name__ == '__main__' :
    #plot where partcle location over time in space
    M31solar = SolarParticles(0,2)
    Sol = M31solar.solarselector()
    #Plot frequency histogram of number of particles at solar radius
    fig, ax = plt.subplots()
    ax.hist(Sol, bins=100, range=(0,10))#,edgecolor='black',color='magenta')
    ax.set(title='Amount of Solarlike particles', xlabel='radius (kpc)', ylabel='frequency') #label graph axes
    plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M31.png')
    
    M31solar = SolarParticles(800,2)
    Sol800 = M31solar.solarselector()
    #Plot frequency histogram of number of particles at solar radius in 12 Gyr
    fig, ax = plt.subplots()
    ax.hist(Sol800, bins=100, range=(0,60),edgecolor='black') #,color='magenta')
    ax.set(title='Solarlike particles', xlabel='radius (kpc)', ylabel='frequency') #label graph axes
    plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M31800.png')
