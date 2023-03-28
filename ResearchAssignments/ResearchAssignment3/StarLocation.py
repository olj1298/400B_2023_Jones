#olivia jones 3.30.23
#code modified from 400B course by G.Besla
"""Fate of Stars at Sun's Location in the M31 Disk.
Using the known values for our Sun's radius and velocity in the Milky Way provided in Martinez-Barbosa+ 2015, 
we can extend this to similar particles that exist in M31. I hope to come to a conclusion on how to select 
sun like particles. I plan to plot circular velocity vs. radius from galactic center."""

import numpy as np
import astropy.units as u
import astropy.constants as const
from ReadFile import Read
from pathlib import Path
from CenterOfMass import CenterOfMass #bring in previous class
import matplotlib.pyplot as plt

#select galactic center of M31
#selecting particles with similar radius as sun in MW (8.2kpc)
#selecting particles with similar velocity as sun in MW (250km/s)
#what snapshots used
#follow particles over time
#when will stop simulation and relevant to which center

#1. Water Maser Distance for the Sun :  R_o = 8.34 kpc   (Reid 2014 ApJ 783) 
#2. GRAVITY Collaboration Distance for the Sun:  R_o = 8.178 kpc   (Abuter+2019 A&A 625)
#3. Value for Distance to Sun listed in Sparke & Gallagher : R_o = 7.9 kpc 
#distances to Galactic Center from the Sun
RoWM = 8.34 * (u.kpc)#kpc
RoG = 9.178 * (u.kpc)#kpc
RoSG = 7.9 * (u.kpc)#kpc
def particalmovement():
    """"Explain.
        Inputs:
            :():
        Returns:
            ::"""
    
    return