import numpy as np 
import astropy.units as u
from ReadFile import Read #use functions from py file
#oliviajones01.26.2023

def ParticleInfo(filename,parttype,partno):
    """Use ReadFile and take inputs like particle type and number
    Inputs:
        :filename(string): address of file to read and access data
        :parttype(integer): value of particle type Dark=1, Disk Star=2, Buldge Star=3
        :partno(float): particle number value
    Returns:
        ::print values for magnitude of distance,velocity, and mass of the particle"""
    
    time,total,data = Read(filename) #assign variables to data in file
    index = np.where(data['type'] == parttype) #find where data in file is particle type user put in
    x = data['x'][index]*u.kpc#data in x column from file. units in kpc 
    y = data['y'][index]*u.kpc #data in y column from file. units in kpc
    z = data['z'][index]*u.kpc #data in z column from file. units in kpc
    vx = data['vx'][index]*u.km/u.s #data in vx column from file. units in km/s   
    vy = data['vy'][index]*u.km/u.s #data in vy column from file. units in km/s
    vz = data['vz'][index]*u.km/u.s #data in vz column from file. units in km/s
    mass = data['m'][index]*u.Msun #data in m column from file. units in M_sun
    distmag = np.sqrt(x[partno-1]**2 + y[partno-1]**2 + z[partno-1]**2) #magnitude of distance of particle
    velmag = np.sqrt(vx[partno-1]**2 + vy[partno-1]**2 + vz[partno-1]**2) #magnitude of velocity of particle
    partmass = mass[partno-1] #mass of particle
    #round to 3 decimals, call variables
    print(f'Magnitude of distance is {np.round(distmag,3)}') #mag of dist specified by user
    print(f'Magnitude of velocity is {np.round(velmag,3)}') #mag of vel specified by user
    print(f'Mass is {partmass}') #mass of particle specified by user
    print(f'Mass is also {np.round(distmag.to(u.lyr),3)}') #convert kpc to lightyr
    return distmag,velmag,partmass