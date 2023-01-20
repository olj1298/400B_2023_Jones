import numpy as np
import astropy.units as u

"""Function to open, read, print, and close file"""

filename = "C:/Users/orang/Downloads/400b/MW_000.txt"
def Read(filename): #filename variable is string of directory address for file to read
    file = open(filename,'r') #open file in directory
    line1 = file.readline() #read line one in file in directory
    label,value = line1.split() #separate values in labels
    time = float(value)*u.Myr  #store time as Myr units
    total = float(value) #store total as number of particles
    file.close() #close file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #dtype=None is line split with white spaces, skip_head=3 is skip 3 first lines, 
    #names=True creates array to store data labels
    #print(data['type'][1])
    return time,total,data #give time, total, data values when asked

Read(filename)