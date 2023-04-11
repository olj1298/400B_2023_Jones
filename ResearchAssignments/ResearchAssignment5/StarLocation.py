#olivia jones 3.30.23
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

#1. Water Maser Distance for the Sun :  R_o = 8.34 kpc (Reid 2014 ApJ 783) 
#2. GRAVITY Collaboration Distance for the Sun:  R_o = 8.178 kpc (Abuter+2019 A&A 625)
#3. Value for Distance to Sun listed in Sparke & Gallagher : R_o = 7.9 kpc 
#distances to Galactic Center from the Sun
RoWM = 8.34 * (u.kpc)#kpc
RoG = 9.178 * (u.kpc)#kpc
RoSG = 7.9 * (u.kpc)#kpc

class SolarParticles: #Class to define COM position and velocity properties of a given galaxy #and simulation snapshot

    def __init__(self,snap,ptype):
        '''Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            Inputs:
            :snap(int): snapshot number of file for particle calculation
            :ptype(integer): 1, 2, or 3. particle type to use for COM calculations'''
        ilbl = '000' + str(snap) #add a string of the filenumber to the value “000”
        ilbl = ilbl[-3:] #remove all but the last 3 digits
        self.filename = 'C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/' + '%s_'%('M31')+ilbl+'.txt'
        self.time, self.total, self.data = Read(self.filename) #Read in the file
        self.index = np.where(self.data['type']==ptype) #Create an index to sort by specified ptype
        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index] #data in m column from file.
        #Create arrays of all position, velocity components
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]    
        
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
        solar = np.array([]) #empty array
        for i in range(len(self.x)): #go through all particles of type specified
            #position and Velocity Vectors
            calcp = np.array([self.x[i], self.y[i], self.z[i]])
            calcv = np.array([self.vx[i], self.vy[i], self.vz[i]])
            #difference of COM to the particles radius and velocity
            radius = np.linalg.norm(calcp - valpvec) #do this as xyz component, move up to compute before
            velocity = np.linalg.norm(calcv - velvec)
            #Add to empty array particles that are needed for this project
            #Using radius from galactic center defined in (Reid 2014 ApJ 783), (Abuter+2019 A&A 625), Sparke & Gallagher
            #using velocity (from Marel+2012)
            if 7.9 < radius < 9.178 and 239 < velocity < 250:
                solar.append(i) #add each fitting particle to empty array
                posintrest = np.array([self.x,self.y,self.z])
                sep = np.array([]) #empty array for radius
                vrel = np.array([]) #empty array for velocity
                for i in solar:
                    r = np.array([self.x[i], self.y[i], self.z[i]]) - valpvec #take difference in radius
                    rv = np.array([self.vx[i],self.vy[i],self.vz[i]]) - velvec #take differnce in velocity
                    sep = np.append(sep, np.linalg.norm(r)) #normalize radius
                    newvrel = np.append(vrel, np.linalg.norm(rv)) #normalize velocity
            #np.where(index) range
            #track particle by spitting out row indices and return r position for same index at snapshot 800    
            #use lab 7 to rotate and get M31 in xyz plane and ring shape so limit on z
        return solar,sep #Return list of all indices with solar properties.
            
if __name__ == '__main__' :
    #plot where partcle location over time in space
    M31solar = SolarParticles(0,2)
    Sol = M31solar.solarselector()
    #Plot frequency histogram of number of particles at solar radius
    fig, ax = plt.subplots()
    ax.hist(Sol, bins=100, range=(0,10))#,edgecolor='black',color='magenta')
    ax.set(title = 'Amount of Solarlike particles', xlabel = 'radius (kpc)', ylabel = 'frequency') #label graph axes
    plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment3/M31.png')
    
    M31solar = SolarParticles(800,2)
    Sol800 = M31solar.solarselector()
    #Plot frequency histogram of number of particles at solar radius in 12 Gyr
    fig, ax = plt.subplots()
    ax.hist(Sol800, bins=100, range=(0,10),edgecolor='black') #,color='magenta')
    ax.set(title = 'Amount of Solarlike particles', xlabel = 'radius (kpc)', ylabel = 'frequency') #label graph axes
    plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment3/M31800.png')

    #ax.semilogy(r,M31H,label = 'Halo Profile')