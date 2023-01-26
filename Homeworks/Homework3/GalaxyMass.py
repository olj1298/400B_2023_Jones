import numpy
import astropy.units as u
from ReadFile import Read #use functions from py file
#oliviajones 2.2.2023

def ComponentMass(filename,parttype):
    """Returns total mass of particle. Store data in a table
    Inputs:
        :filename(string): address of file to read and access data
        :parttype(integer): value of particle type Dark=1, Disk Star=2, Buldge Star=3     
    Returns: 
        ::Mass in units of 10^12 M_sun"""


#directory for file to read
filename = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework2/MW_000.txt"
ComponentMass(filename,2)    