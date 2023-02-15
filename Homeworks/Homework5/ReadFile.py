import numpy as np
import astropy.units as u
from pathlib import Path
#oliviajones01.26.2023

def Read(filename): 
    """Function to open, read, print, and close file
    Inputs:
        :filename(string): address of file to read and access data
    Returns:
        ::time, total, data values"""
    file = open(filename,'r') #open file in directory
    line1 = file.readline() #read line one in file in directory
    label,value = line1.split() #separate values in labels
    time = float(value)*u.Myr  #store time as Myr units
    
    line2 = file.readline()
    label2,value2 = line2.split()
    total = float(value2) #store total as number of particles
    file.close() #close file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #dtype=None is line split with white spaces, skip_head=3 is skip 3 first lines, 
    #names=True creates array to store data labels
    #print(data['type'][1])
    return time,total,data #give time, total, data values when asked

#Read('MW_000.txt')

def main():
#define the file path to be used to open the file
    filename=Path('C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework5/MW_000.txt')
    
    #call Read function to produce a numpy array of the data within the file
    Read(filename)
    
main()