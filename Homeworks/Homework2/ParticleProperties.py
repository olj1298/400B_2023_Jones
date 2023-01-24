import numpy as np 
import astropy.units as u
from ReadFile import Read #use functions from 
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
    parttype = data['type'] 
    x = data['x']*u.kpc #data in x column from file. units in kpc 
    y = data['y']*u.kpc #data in y column from file. units in kpc
    z = data['z']*u.kpc #data in z column from file. units in kpc
    vx = data['vx']*u.kms #data in vx column from file. units in km/s   
    vy = data['vy']*u.kms #data in vy column from file. units in km/s
    vz = data['vz']*u.kms #data in vz column from file. units in km/s
    partno = data['i']
    index = np.where(data['x']>2)
    xnew = data['x'][index]*u.kpc
    
    #round to 3 units
    print(f'Magnitude of distance is {data.distance} kpc')
    print(f'Magnitude of velocity is {data.speed} km/s')
    print(f'Mass is {data.mass} M_Sun')
    round(x,3)