import numpy as np
import pandas as pd
import astropy.units as u
from ReadFile import Read #use functions from py file
from tabulate import tabulate #table functions
#oliviajones 2.2.2023

def ComponentMass(filename,parttype):
    """Returns total mass of particle type. Store data in a table
    Inputs:
        :filename(string): address of file to read and access data
        :parttype(integer): value of particle type Halo=1, Disk=2, Buldge=3     
    Returns: 
        ::Mass in units of 10^12 M_sun"""
    
    time,total,data = Read(filename) #assign variables to data in file
    index = np.where(data['type'] == parttype) #find where data in file is particle type user put in
    mass = data['m'][index] #data in m column from file.
    masssum = (np.sum(mass)*10e10*u.Msun) #sum all particle type masses in m column units in 10e12 M_sun
    return masssum/100 #totalmass

#directory for file to read
MWfile = "C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/MW_000.txt"
M31file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M31_000.txt"
M33file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/ResearchAssignments/ResearchAssignment5/M33_000.txt"

totmassMW1 = ComponentMass(MWfile,1) #Mass of all Halo particles in Milky Way
#print(f'Mass of halo particles in Milky Way is {totmassMW1:.2e}') #round answer to two places
totmassMW2 = ComponentMass(MWfile,2) #Mass of all disk particles in Milky Way   
#print(f'Mass of disk particles in Milky Way is {totmassMW2:.2e}') #round answer to two places
totmassMW3 = ComponentMass(MWfile,3) #Mass of all Buldge particlces in Milky Way
#print(f'Mass of buldge particles in Milky Way is {totmassMW3:.2e}') #round answer to two places
totmassMW = totmassMW1+totmassMW2+totmassMW3 #sum of all particle types
#print(f'Total mass of the Milky Way is {totmassMW:.2e}') #round answer to two places
fbarMW = (totmassMW2+totmassMW3)/totmassMW #fbar is ((disk + buldge) particles mass/totalmass of all particles)
#print(f'f_bar of M33 is {fbarMW:.3E}') #round answer to two places

totmassM311 = ComponentMass(M31file,1)#Mass of all Halo particles in M31
#print(f'Mass of halo particles in M31 is {totmassM311:.2e}') #round answer to two places
totmassM312 = ComponentMass(M31file,2)#Mass of all Disk particles in M31
#print(f'Mass of disk particles in M31 is {totmassM312:.2e}') #round answer to two places
totmassM313 = ComponentMass(M31file,3)#Mass of all Buldge particlces in M31
#print(f'Mass of buldge particles in M31 is {totmassM313:.2e}') #round answer to two places
totmassM31 = totmassM311+totmassM312+totmassM313 #sum of all particle types
#print(f'Total mass of the M31 is {totmassM31:.3E}') #round answer to two places
fbarM31 = (totmassM312+totmassM313)/totmassM31 #fbar is ((disk + buldge) particles mass/totalmass of all particles)
#print(f'f_bar of M33 is {fbarM31:.3E}') #round answer to two places

totmassM331 = ComponentMass(M33file,1)#Mass of all Halo particles in M33
#print(f'Mass of halo particles in M33 is {totmassM331:.2e}') #round answer to two places
totmassM332 = ComponentMass(M33file,2)#Mass of all Disk particles in M33  
#print(f'Mass of disk particles in M33 is {totmassM332:.2e}') #round answer to two places
#print(f'Mass of buldge particles in M33 is zero since no buldge particles.')
totmassM33 = totmassM331+totmassM332 #sum of all particle types
#print(f'Total mass of the M33 is {totmassM33:.2e}') #round answer to three places
fbarM33 = totmassM332/totmassM33 #fbar is (disk particles mass/totalmass of all particles)
#print(f'f_bar of M33 is {fbarM33:.3E}') #round answer to two places

localmass = totmassMW+totmassM31+totmassM33 #sum of particles in each galaxy in group
#print(f'Mass of Local Group is {localmass:.2e}') #round answer to two places

#assign data
mydata =[["Milky Way", np.round(totmassMW1/10e12,3), np.round(totmassMW2/10e12,3), np.round(totmassMW3/10e12,3), 
          np.round(totmassMW/10e12,3), np.round(fbarMW,3)], #row 1 
            ["M31", np.round(totmassM311/10e12,3), np.round(totmassM312/10e12,3), np.round(totmassM313/10e12,3), 
             np.round(totmassM31/10e12,3), np.round(fbarM31,3)], #row 2
            ["M33", np.round(totmassM331/10e12,3), np.round(totmassM332/10e12,3), "zero", 
             np.round(totmassM33/10e12,3), np.round(fbarM33,3)], #row 3
            ["Local Group", "zero", "zero", "zero", np.round(localmass/10e12,3), "zero"] #row 4
        ]
#create header
head = ["Galaxy Name", "Halo Mass (10e12 M_sun)", "Disk Mass (10e12 M_sun)", 
        "Buldge Mass (10e12 M_sun)", "Total (10e12 M_sun)", "fbar"]

#display table
#print(tabulate(mydata, headers=head, tablefmt="grid"))