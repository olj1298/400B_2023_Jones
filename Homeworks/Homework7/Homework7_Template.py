#Homework 7 Template
#Rixin Li & G . Besla
#Make edits where instructed - look for "****", which indicates where you need to add code.
#oliviajones3.23.23 

#import necessary modules
#numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
#matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
#astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
#import Latex module so we can display the results with symbols
from IPython.display import Latex
#**** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass #bring in previous class
#**** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass
#M33AnalyticOrbit
from OrbitCOM import VectorDifference
import os

class M33AnalyticOrbit:
    """Calculate the analytical orbit of M33 around M31"""
    
    def __init__(self,filename):
        """Velocity and position calculation.
        Inputs:
            :filename(string): directory address to store integrated orbit"""

        #gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        #store the output file name
        self.filename=filename
        #get the current pos/vel of M33 
        #create an instance of the  CenterOfMass class for M33 
        COM_M33=CenterOfMass('M33_000.txt',2)
        #store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_delta=.1
        M33_voldec=4
        COM_pos_M33=COM_M33.COM_P(M33_delta,M33_voldec)
        
        #store the velocity VECTOR of the M33 COM (.value to get rid of units)
        COM_vel_M33=COM_M33.COM_V(COM_pos_M33[0],COM_pos_M33[1],COM_pos_M33[2]).value
        
        #get the current pos/vel of M31 
        #**** create an instance of the  CenterOfMass class for M31 
        COM_M31=CenterOfMass('M31_000.txt',2)
        #**** store the position VECTOR of the M31 COM (.value to get rid of units)
        delta=.1
        voldec=2
        COM_pos_M31=COM_M31.COM_P(delta,voldec)
        #**** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        COM_vel_M31=COM_M31.COM_V(COM_pos_M31[0],COM_pos_M31[1],COM_pos_M31[2]).value
        
        #store the DIFFERENCE between the vectors posM33 - posM31
        #*** create two VECTORs self.r0 and self.v0 and have them be the
        #relative position and velocity VECTORS of M33
        self.r0=(COM_pos_M33-COM_pos_M31).value
        self.v0=(COM_vel_M33-COM_vel_M31)
        
        #get the mass of each component in M31
        
        #disk
        #**** self.rdisk = scale length (no units)
        self.rdisk = 5 #kpc
        #**** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk=(ComponentMass('M31_000.txt',2)*1e12).value #from HW3 
        
        #bulge
        #**** self.rbulge = set scale length (no units)
        self.rbulge = 1 #kpc
        #**** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = (ComponentMass('M31_000.txt',3)*1e12).value #from HW3
        
        #halo
        #**** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 61 #kpc
        #**** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = (ComponentMass('M31_000.txt',1)*1e12).value #from HW3
    
    def HernquistAccel(self,M,r_a,r): # it is easiest if you take as an input the position VECTOR 
        """Calculate acceleration vector using the Hernquist model for the halo and bulge masses.
        Inputs:
            :self():
            :M(float): total halo or bulge mass of M31
            :r_a(float): scale length of M31 for either the bulge or halo
            :r(array): relative COM position vector of M31 and M33
        Returns:
            ::This function returns the acceleration vector from a Hernquist potential"""
        
        #Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        
        #Store the Acceleration
        Hern = -self.G*M/(rmag*(r_a+rmag)**2)*r #follow the formula in the HW instructions
        #NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        #use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    def MiyamotoNagaiAccel(self,M,r_d,r): #it is easiest if you take as an input a position VECTOR  r 
        """Acceleration vector using the Miyamoto-Nagai profile based on the disk properties.
        Inputs:
            :self():
            :M(float): total halo or bulge mass of M31
            :r_d(float): scale length for the disk of M31
            :r(array): relative COM position vector of M31 and M33
        Returns:
            ::This function returns the acceleration vector from a Miyamoto-Nagai profile"""

        z_d = self.rdisk/5
        R = np.sqrt(r[0]**2+r[1]**2)
        B = r_d+np.sqrt(r[2]**2+z_d**2)
        #define the constants in all accelerations
        numerator = -self.G*M
        denominator = (R**2+B**2)**1.5
        
        ax = numerator*r[0]/denominator
        ay = numerator*r[1]/denominator
        az = numerator*r[2]*B/(denominator*np.sqrt(r[2]**2+(r_d/5)**2))
        
        accel = np.array([ax,ay,az])
        
        return accel
        
        #Acceleration **** follow the formula in the HW instructions
        #AGAIN note that we want a VECTOR to be returned (see Hernquist instructions)
        #this can be tricky given that the z component is different than in the x or y directions. 
        #we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        #differences in the z direction:
        #multiply the whle thing by : np.array([1,1,ZSTUFF]) 
        #where ZSTUFF are the terms associated with the z direction
        #the np.array allows for a different value for the z component of the acceleration
     
    def M31Accel(self,COM_M31): #input should include the position vector, r
        """Total acceleration vector from the sum of the component vectors.
            Inputs:
                :self():
                :COM_M31(array): relative COM position vector of M31 and M33
            Returns:
                ::3D vector of the total acceleration"""

        #Call the previous functions for the halo, bulge and disk
        #these functions will take as inputs variable we defined in the initialization of the class like self.rdisk etc. 
        a_halo=self.HernquistAccel(self.Mhalo,self.rhalo,COM_M31)
        a_bulge=self.HernquistAccel(self.Mbulge,self.rbulge,COM_M31)
        a_disk=self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,COM_M31)
        
        sum_a=a_halo+a_bulge+a_disk    
        
        return sum_a #.value #return the SUM of the output of the acceleration functions - this will return a VECTOR 
    
    def LeapFrog(self,dt,r,v): #take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """change in COM position and velocity of M33 by taking a step forward (or backward if dt<0) in time.
            Inputs:
                :self():
                :dt(float): time interval for integration
                :r(array): COM xyz position of M33 relative to M31
                :v(array): xyz velocity of the COM of M33 relative to M31
            Returns:
                ::COM position and velocity vector of M33 after taking a step dt in time"""
        
        #predict the position at the next half timestep
        rhalf = r+v*dt/2
        
        #predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v+self.M31Accel(rhalf)*dt
        
        #predict the final position using the average of the current velocity and the final velocity
        #this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        #next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf+vnew*dt/2
        
        return rnew,vnew #return the new position and velcoity vectors
    
    def OrbitIntegration(self, t0, dt, tmax):
        """orbit of M33 by integrating forward (or backward if dt<0) in time
            Inputs:
                :self():
                :t0(float): starting time
                :dt(float): time interval
                :tmax(float): final time"""

        #initialize the time to the input starting time
        t = t0
        
        #initialize an empty array of size: rows int(tmax/dt)+2, columns 7
        orbit = np.zeros((int(tmax/dt)+2,7))
        
        #initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        #this above is equivalent to orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        #initialize a counter for the orbit.  
        i = 1 #since we already set the 0th values, we start the counter at 1
        r,v = self.r0,self.v0
        
        #start the integration (advancing in time steps and computing LeapFrog at each step)
        while t<tmax: #as long as t has not exceeded the maximal time 
            
            #advance the time by one timestep, dt
            t+=dt
            
            #store the new time in the first column of the ith row
            orbit[i,0] = t
            
            #advance the position and velocity using the LeapFrog scheme
            #remember that LeapFrog returns a position vector and a velocity vector  
            #as an example, if a function returns three vectors you would call the function and store 
            #the variable like:a,b,c = function(input)
            rnew,vnew=self.LeapFrog(dt,r,v)
         
            #store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            #TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : A[n, 5:8] 
            #where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i,1:4] = rnew
            orbit[i,4:8] = vnew
            
            #store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit            
            #update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i+=1
            r = rnew
            v = vnew
            
        #write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

if __name__ == '__main__' :
    #set the font size for all figures
    plt.rcParams.update({'font.size': 18})
    #set the filename for thr file that will hold the M33 positions
    filename='M33_integrated_pos.txt'
    #initialize class that will do the orbital integration
    number=M33AnalyticOrbit(filename).OrbitIntegration(t0=0,dt=0.1,tmax=10)
    #load M33 integration data into array
    data=np.genfromtxt(os.path.abspath(filename),dtype=None,names=True,skip_header=0)
    #create M33 distance, velocity, and time arrays
    r=np.sqrt(data['x']**2+data['y']**2+data['z']**2)
    v=np.sqrt(data['vx']**2+data['vy']**2+data['vz']**2)
    t=data['t']
    
    '''import data from galaxy merger'''
    #open M31 data fron the galaxy merger simulation
    M31_file=file=open('Orbit_M31_0_800.txt','r')
    M31_data=np.genfromtxt(os.path.abspath(M31_file),dtype=None,names=True,skip_header=0)
    #open M33 data fron the galaxy merger simulation
    M33_file=file=open('Orbit_M33_0_800.txt','r')
    M33_data=np.genfromtxt(os.path.abspath(M33_file),dtype=None,names=True,skip_header=0)
    #find the relative position and velocity of M33 and M31 from the merger simulation
    merge_pos=VectorDifference([M33_data['x'],M33_data['y'],M33_data['z']],[M31_data['x'],M31_data['y'],M31_data['z']])
    merge_vel=VectorDifference([M33_data['vx'],M33_data['vy'],M33_data['vz']],[M31_data['vx'],M31_data['vy'],M31_data['vz']])
    merge_t=M33_data['t']
    
    #plot
    fig,ax=plt.subplots(1,2,figsize=(20,10))
    #plot M33 and M31 relative position as a function of time for the integration
    #done by this script and for the merger simulation
    ax[0].set(title='M33 and M31 Relative Distance', xlabel='Time (Gyr)', ylabel='Distance (kpc)')
    ax[0].plot(t[:len(t)-1],r[:len(r)-1],label='Integration')
    ax[0].plot(merge_t,merge_pos,label='Simulation')
    ax[0].legend()
    #plot M33 and M31 relative velocity as a function of time for the integration
    #done by this script and for the merger simulation  
    ax[1].set(title='M33 and M31 Relative Velocity', xlabel='Time (Gyr)', ylabel='Velocity (km/s)')
    ax[1].plot(t[:len(t)-1],v[:len(v)-1],label='Integration')
    ax[1].plot(merge_t,merge_vel,label='Simulation')
    ax[1].legend()
    plt.savefig('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework7/M33distance_velocity.png')

''''Answer HW Questions'''
    
print(f"Question 2) The plots diverge after about 1.5 Gyr. The analytic orbit is one that is circular but digresses to a sin curve over time.,\
      since we don't assume change to the system, the original orbit continues.")
    
print(f"Question 3) Interactions with the Milky Way from merger is not included in the simulation. \
        This assumes M31 is also constant over time. There is assumption that galaxies are point masses, not interaction between particles.\n")
          
print(f"Question 4) To include the Milky Way, the forces and acceleration \
        between all 3 galaxies would need to be computed, making the computation more complex.")