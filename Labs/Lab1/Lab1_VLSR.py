#In Class Lab 1
#Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023
import numpy as np 
import astropy.units as u 
from astropy import constants as const #import astropy constants

#Part A:  The Local Standard of Rest
#Proper motion of Sgr A* from Reid & Brunthaler 2004 $\mu = 6.379$ mas/yr 
#Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
#$v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$

#a) Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
#The function should take as input: the solar radius (R$_o$), the proper motion (mu)
#and the peculiar motion of the sun in the $v_\odot$ direction.
mu = 6.379 #proper motion of sun
vsun=12.24* (u.km/u.s) #velocity of sun in v direction

#Compute V$_{LSR}$ using three different values R$_o$: 
#1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
#2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
#3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
#distances to Galactic Center from the Sun
RoWM = 8.34 * (u.kpc)#kpc
RoG = 9.178 * (u.kpc)#kpc
RoSG = 7.9 * (u.kpc)#kpc

def VLSR(r0, mu = 6.379, vsun=12.24* (u.km/u.s)):
    """Compute velocity at local standard of rest. VLSR = 4.74*mu*r0 - vsun
    Inputs:
        :r0(astropy quantity): distance from sun to galactic center in kpc
        :mu(float): proper motion of Sag A* in mas/yr. Default from Reid & Brunthaler 2004
        :vsun(astropy quantity): peculiar motion of the sun. Default is from Schonrich + 2010
    Returns:
        ::VLSR in km/s"""
    return 4.74 * (mu) * (r0/u.kpc) * (u.km/u.s) - vsun

VLSRreid = VLSR(RoWM,mu,vsun) #VLSR using Reid
print(f'VLSR using Reid+2014 is {VLSRreid}')

VLSRgravity = VLSR(RoG,mu,vsun) #VLSR using GRAVITY
print(f'VLSR using GRAVITY collab is{VLSRgravity}')

VLSRsg = VLSR(RoSG,mu,vsun) #VLSR using Sparke and Gallagher
print(f'VLSR using Sparke+ is{VLSRsg}')

#b)compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
#Note that 1 km/s $\sim$ 1kpc/Gyr

def TorbSun(R, V):
    """Computes orbital period of the Sun. T = 2*pi*R / V
    Inputs:
        :R(astropy quantity): distance to galactic center in kpc
        :V(astropy quantity): velocity in km/s of sun in v direction
    Returns:
        ::astropy quantity of orbital period in Gyr"""
    VkpcGyr = V.to(u.kpc/u.Gyr)
    T = 2*np.pi*R/VkpcGyr    
    return T

VsunPeculiar = 12.23*u.km/u.s #peculiar motion
VSun = VLSRgravity +  VsunPeculiar #find vel of Sun in v direction using VLSR and V peculiar

T_Grav = TorbSun(RoG,VSun)

#c)Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)
#Age of Universe / Orbital Period
Age = 13.8*u.Gyr
print(f'Number of rotations about GC over life of universe is {Age/T_Grav} rotations') 

#Part B Dark Matter Density Profiles
#a)Try out Fitting Rotation Curves [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
#Point Mass = 0, Star Mass = 80, Dark Matter = 12

#b)In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
#Recall that for the Isothermal sphere :
#$\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
#Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
#What about at 260 kpc (in units of  M$_\odot$) ?

#density profile rho = VLSR**2 / (4*pi*G*R**2)
#mass = integrate rho dV
#     = rho 4*pi*r**2 dr
#     = VLSR**2/G/(4*pi*r**2) * (4*pi*r**2) dr
#     = VLSR**2 / G*r

Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun) #Gravitational Constant

def MassIso(r,VLSR):
    """Compute dark matter enclosed in distance assuming isothermal sphere model for dark matter.
    Inputs:
        :r(astropy quantity): distance to Galactic Center in kpc
        :VLSR(astropy quantity): velocity of the Local Standard of Rest in km/s"""

    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) #convert km/s to kpc/Gyr
    M = VLSRkpcGyr ** 2 / Grav * r #mass for isothermal sphere
    return M

MIsoSolar = MassIso(RoG,VLSRgravity) #Mass using GRAVITY data
print(MIsoSolar)

print(f'Mass using GRAVITY data {MIsoSolar:.2e}')

#mass in 260kpc
MIso260 = MassIso(260*u.kpc, VLSRgravity)
print(f'Mass inside 260 kpc is {MIso260:.2e}')

#c)The Leo I satellite is one of the fastest moving satellite galaxies we know. 
#It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
#If we assume that Leo I is moving at the escape speed: $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
#and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the min. mass of the Milky Way (in units of M$_\odot$) ?  
#How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

Vtot = 196 *u.km/u.s

def Leosprint():
    return