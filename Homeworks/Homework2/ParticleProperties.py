import numpy as np 
import astropy 
import astropy.units as u
from ReadFile import Read #use functions from 

"""Use ReadFile and take inputs like particle type and number"""

def ParticleInfo(filename,parttype,partno):
    print(f"Magnitude of distance is {Read.distance} kpc")
    print(f"Magnitude of velocity is {rf.speed} km/s")
    print(f"Mass is {rf.mass} M_Sun")
    