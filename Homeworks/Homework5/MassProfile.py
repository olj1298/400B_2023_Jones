#olivia jones 2.16.23
#import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read #bring in readfile.py
from CenterOfMass import CenterOfMass #bring in previous class
import matplotlib.pyplot as plt

#In this assignment you will determine the mass distribution of each galaxy at SnapNumber and use this to determine each galaxy’s rotation curve.
class MassProfile:

    def __init__(self, galaxy, snap):
        '''In this assignment you will determine the mass distribution of each galaxy at SnapNumber
        and use this to determine each galaxy's rotation curve. 
            Inputs:
            :self(function within class):
            :galaxy(string): Galaxy Name such as MW, M31 and M33
            :snap(integer): Snapshot number 1, 2, etc.'''
        
        ilbl = '000' + str(snap) #add a string of the filenumber to the value “000”
        #remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename='%s_'%(galaxy)+ilbl+'.txt'
        #read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                      
        #store the mass, positions, velocities, galaxy name from file
        self.m = self.data['m'] #data in m column from file. 
        self.x = self.data['x']*u.kpc #data in x column from file.  
        self.y = self.data['y']*u.kpc #data in y column from file.  
        self.z = self.data['z']*u.kpc #data in z column from file.  
        self.gname = galaxy #galaxy name
    
    def MassEnclosed(self, ptype,rarray):
        """Computes the mass enclosed within a given radius of the COM position for a specified galaxy and a 
        specified component of that galaxy. Calculate the enclosed masses for an array of radii and returns an array
        of masses in Msun. 
            Inputs:
                :self(function within class):
                :ptype(): particle type
                :rarray(): array of radii as mag., not vector
            Returns:
                ::mass enclosed(array): units in Msun"""
                
        GalCOM = CenterOfMass(self.filename, ptype) #Create a Center of mass object for COM Postion
        GalCOMp = GalCOM.COM_P(0.1) #Call COM_P
        indpart = np.where(self.data["type"] == ptype) #For particle type that exists in galaxy
        #define the x, y, and z center of masses
        x_COM = GalCOMp[0] #0=x
        y_COM = GalCOMp[1] #1=y
        z_COM = GalCOMp[2] #2=z
        #positions of the particles rel. to the COM
        x_new = self.x[indpart]-x_COM
        y_new = self.y[indpart]-y_COM
        z_new = self.z[indpart]-z_COM
        m_new = self.m[indpart]
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) #radius of each particle from the center of the galaxy
        mlist=[] #empty list for mass values

        for radius in rarray: #for radius in array of radii
            rindex = np.where(r_new < radius) #get particles inside radius
            minr = m_new[rindex] #get mass of particles in radius
            mlist.append(sum(minr)) #sum mass of all particles in radius
         
        massarray = np.array(mlist) #change from list to array
        
        return massarray*u.M_sun*1e10 #add units and magnitude
    
    def MassTotal(self,rarray):
        """Calls MassEnclosed to compute the mass enclosed within the radius array for each
        particle type (bulge, disk and halo).
            Inputs:
                :self(function within class):
                :rarray(array): 1D array of radii
            Returns:
                ::array of masses, enclosed mass (bulge+disk+halo) at each radius of the input array.units in Msun"""

        Hmass = self.MassEnclosed(1,rarray) #call mass array enclosed in halo
        Dmass = self.MassEnclosed(2,rarray) #call mass array enclosed in disk
        
        if self.gname != 'M33': #exception for two galaxys with bulge mass
            Bmass = self.MassEnclosed(3,rarray) #include bulge mass
            Tmass = Hmass+Dmass+Bmass #add halo,disk,bulge for total mass
        else: #exception for galaxy without bulge mass
            Tmass = Hmass+Dmass #only add halo and disk mass
        
        return Tmass #units in Msun
    
    def HerquistMassProfile(self, radius, a, Mhalo):
        """Compute the mass enclosed within a given radius using the theoretical Herquist profile.
                p(r) = (M*a / 2 pi r) * (1/(r+a)**3) 
                M(r) = (Mhalo* r**2) / (a+r)**2
            Inputs:
                :self(function within class):
                :radius(integer): radius of galaxy in kpc
                :a(integer): scale factor
                :Mhalo(integer): Halo mass in Msun
            Returns:
                ::Halo mass in units of Msun"""

        M = (Mhalo * radius**2) / (a+radius)**2 #calc mass
        p = (M * a / (2 * np.pi * radius)) * (1/(radius+a)**3) #hernquist profile
        
        Hmassprint = print(f"Halo Mass is: {np.round(p,2)*u.Msun}")
        return Hmassprint
    
    def CircularVelocity(self,ptype,rarray):
        """The circular speed is computed using the Mass enclosed at each radius, assuming spherical symmetry.
                G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
            Inputs:
                :self(function within class):
                :ptype(integer): particle type
                :rarray(array): 1D array of radii
            Returns:
                ::array of circular speeds in units of km/s, rounded to two decimal places."""
        
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #Grav const. conversion
        marray = self.MassEnclosed(ptype,rarray) #Call mass at radius for particles
        vel = np.sqrt(G*marray/rarray) #Orbital velocity equation
        
        Circprint = print(f"velocity is: {np.round(vel,2)*u.km/u.s}") #Round and get in units for vel.
        return Circprint

    def CircularVelocityTotal(self,rarray):
        """The total circular velocity is NOT just the circular velocity of each individual galaxy component summed together.
            Inputs:
                :self(function within class):
                :rarray(array): 1D array of radii
            Returns:
                ::array of circular velocity (in units of km/s) representing the total Vcirc
                created by all the galaxy components (bulge+disk+halo) at each radius of the input
                array."""

        Tmarray = self.MassEnclosedTotal(rarray) #Call total mass of particles per radius
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #Grav const. conversion
        Tvcirc = np.sqrt(G*Tmarray/rarray) #Orbital velocity equation
        
        Cvelprint = print(f"total velocity: {np.round(Tvcirc,2)*u.km/u.s}") #round and put in units for vel.
        return Cvelprint
    
    def HerquistVCirc(self,radius,a,Mhalo):
        """Computes the circular speed using the Hernquist mass profile.
            Inputs:
                :self(function within class):
            Returns:
                ::circular speed in units of km/s, rounded to two decimal places."""
          
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #Grav const. conversion
        Hernquist = self.HernquistMass(radius,a,Mhalo) #Call Herquist profile with mass and radius
        VCirc = np.sqrt(G*Hernquist*u.Msun/(radius)) #Orbital velocity equation with Herquist profile
        
        HCircprint = print(f"Herquist VCirc: {np.round(VCirc,2)*u.km/u.s}") #round and put in units for vel.                 
        return HCircprint
    
if __name__ == '__main__' :
    #Section 8
    #1. Create a plot of the mass profile (mass enclosed as a function of radius) of each component of the MW to a radius of 30 kpc.
    #2. You will need to define an array of radii. Don’t start the first element at 0 (say 0.1 instead).
    #3. Use different line types and/or colors to indicate different components. If you get stuck with the plotting look at
    #the solutions for InClassLab1 on the class GitHub repo.
    #4. Y-axis should be in log (plt.semilogy).
    #5. Overplot the sum of all components (MassEnclosedTotal)
    #6. Determine the best fitting Hernquist profile:
        #Use the total dark matter mass of each galaxy as Mhalo.
        #You will need to find the scale radius (a) where the Hernquist profile best matches the
        #plotted mass profile as a function of radius. Do this by plotting the Hernquist profile
        #(starting with a guess for a) on top of the simulated dark matter halo mass profile.
        #Then keep changing a until the mass distributions reasonably agree. Alternatively, you
        #could print out values for the simulated enclosed mass at 30 kpc and change a until
        #you get a reasonable agreement with the Hernquist profile.
    #7. Use a different line-type or color for the best fit Hernquist profile and include text that
    #states the best fit scale length.
    
    #directory address for Galaxy files
    #MWfile = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework5/MW_000.txt"
    #M31file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework5/M31_000.txt"
    #M33file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework5/M33_000.txt"
    
    r = np.arange(0.25,30.5,0.5)*u.kpc #radius to get mass and vel.
    
    """MW"""
    MW = MassProfile('MW',0) #COM object for milky way
    MWH = MW.MassEnclosed(1,r) #mass profile of MW Halo
    MWD = MW.MassEnclosed(2,r) #mass profile of MW Disk
    MWB = MW.MassEnclosed(3,r) #mass profile MW Bulge
    MWT = MW.MassEnclosedTotal(r) #MW Mass Total
    MassMWH = max(MWH)/u.Msun #Mass of MW halo
    Hern = [] #empty array for Herquist profile masses
    scaleMW = 18 #scale factor of MW
    
    for radius in r: #find the mass according to hernquist profile at each radii
        MWHern = MW.HernquistMass(radius/u.kpc,scaleMW,MassMWH)
        Hern.append(MWHern)
    
    #plot MW mass profiles on log graphs
    fig,ax = plt.subplots(figsize = (10,10))
    ax.semilogy(r,MWH,label = 'Halo Profile')
    ax.semilogy(r,MWD,label = 'Disk Profile')
    ax.semilogy(r,MWB,label = 'Bulge Profile')
    ax.semilogy(r,MWT,label = 'Total Profile')
    ax.semilogy(r,Hern,linestyle = 'dashed',label = 'Hernquist Profile')
    legend = ax.legend() #create legend
    ax.set(title = 'MW Mass Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Mass Enclosed ($M_{\odot}$))') #label graph axes
    print(f"The best scale length for the MW is: {scaleMW*u.kpc}")
    
    #8. Repeat for M33 and M31.
    
    """M31"""
    M31 = MassProfile('M31',0) #COM object for M31
    M31H = M31.MassEnclosed(1,r) #mass profile of M31 Halo
    M31D = M31.MassEnclosed(2,r) #mass profile of M31 Disk
    M31B = M31.MassEnclosed(3,r) #mass profile of M31 Bulge
    M31T = MW.MassEnclosedTotal(r) #M31 mass total
    MassM31H = max(M31H)/u.Msun #mass of M31 halo
    Hern=[] #empty array for Herquist profile masses
    scaleM31=15 #scale factor of M31
    
    for radius in r: #mass according to hernquist profile at each radii
        M31Hern=M31.HernquistMass(radius/u.kpc,scaleM31,MassM31H)
        Hern.append(M31Hern)
    
    #plot M31 mass profiles on log graphs
    fig,ax=plt.subplots(figsize = (10,10))
    ax.semilogy(r,M31H,label = 'Halo Profile')
    ax.semilogy(r,M31D,label = 'Disk Profile')
    ax.semilogy(r,M31B,label = 'Bulge Profile')
    ax.semilogy(r,M31T,label = 'Total Profile')
    ax.semilogy(r,Hern,linestyle = 'dashed',label = 'Hernquist Profile')
    legend=ax.legend() #create legend
    ax.set(title = 'M31 Mass Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Mass Enclosed ($M_{\odot}$))') #label graph axes
    print(f"The best scale length for M31 is: {scaleM31*u.kpc}")
    
    """M33"""
    M33 = MassProfile('M33',0) #COM object for M33
    M33H = M33.MassEnclosed(1,r) #mass profile of M33 Halo
    M33D = M33.MassEnclosed(2,r) #mass profile of M33 Disk
    M33T = MW.MassEnclosedTotal(r) #M33 mass total
    MassM33H = max(M33H)/u.Msun #mass of M33 halo
    Hern = [] #empty array for Herquist profile masses
    scaleM33 = 10 #specify scale factor of M33
    
    for radius in r: #mass according to hernquist profile at each radiii
        M33Hern = M33.HernquistMass(radius/u.kpc,scaleM33,MassM33H)
        Hern.append(M33Hern)
    
    #plot M33 mass profiles on log graphs
    fig,ax = plt.subplots(figsize = (10,10))
    ax.semilogy(r,M33H,label = 'Halo Profile')
    ax.semilogy(r,M33D,label = 'Disk Profile')
    ax.semilogy(r,M33T,label = 'Total Profile')
    ax.semilogy(r,Hern,linestyle = 'dashed',label = 'Hernquist Profile',color = 'red')
    legend=ax.legend() #create legend
    ax.set(title='M33 Mass Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Mass Enclosed ($M_{\odot}$))')#label graph axes
    print(f"The best scale length for M33 is: {scaleM33*u.kpc}")
    
    #9 Plot the Rotation Curve for each Galaxy
    
    #Section 9
    #1. Create a plot of the Rotation Curve of each galaxy - i.e. the circular speed of each
    #galaxy component to a radius of 30 kpc.

    #2. You will need to define an array of radii. Don’t start the first element at 0 (say 0.1
    #instead).

    #3. Use different line types and/or colors to indicate different components.

    #4. Overplot the total circular speed for all galaxy component (CircularVelocityTotal).

    #5. Overplot the best fit theoretical Hernquist circular speed, using the scale radius you
    #determined from the mass enclosed plot.
    
    """MW Rotation Curves"""
    MWHv = MW.CircularVelocity(1,r) #vel. profile of MW Halo
    MWDv = MW.CircularVelocity(2,r) #vel. profile of MW Disk
    MWBv = MW.CircularVelocity(3,r) #vel. profile of MW Bulge
    MWTv = MW.CircularVelocityTotal(r) #Total vel. of MW
    Hernv = [] #empty array for Herquist vel. profile
    for radius in r: #vel. according to hernquist profile at each radius
        MWHernv = MW.HernquistVCirc(radius,scaleMW*u.kpc,MassMWH)
        Hernv.append(MWHernv*u.s/u.km)
    
    #plot MW velocity profiles on log graphs
    fig,ax=plt.subplots(figsize = (10,10))
    ax.scatter(r,MWHv,label = 'Halo Profile')
    ax.scatter(r,MWDv,label = 'Disk Profile')
    ax.scatter(r,MWBv,label = 'Bulge Profile')
    ax.scatter(r,MWTv,label = 'Total Profile')
    ax.plot(r,Hernv,linestyle = 'dashed',label = 'Hernquist Profile',color = 'red')
    legend=ax.legend(loc = 8) #create legend
    ax.set(title='MW Velocity Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Velocity ($ms^{-1}$))') #label graph axes

    #6. Do the same for M31 and M33.
    
    """M31 Rotation Curve"""
    M31Hv = M31.CircularVelocity(1,r) #vel. profile M31 Halo
    M31Dv = M31.CircularVelocity(2,r) #vel. profile M31 Disk
    M31Bv = M31.CircularVelocity(3,r) #vel. profile M31 Bulge
    M31Tv = M31.CircularVelocityTotal(r) #Total vel. profile for M31
    Hernv = [] #empty array for Herquist vel. profile
    for radius in r: #vel. according to hernquist profile at each radius
        M31Hern = M31.HernquistVCirc(radius,scaleM31*u.kpc,MassM31H)
        Hernv.append(M31Hern*u.s/u.km)
    
    #plot M31 velocity profiles on log graphs
    fig,ax = plt.subplots(figsize = (10,10))
    ax.scatter(r,M31Hv,label = 'Halo Profile')
    ax.scatter(r,M31Dv,label = 'Disk Profile')
    ax.scatter(r,M31Bv,label = 'Bulge Profile')
    ax.scatter(r,M31Tv,label = 'Total Profile')
    ax.plot(r,Hernv,linestyle = 'dashed',label = 'Hernquist Profile',color = 'red')
    legend = ax.legend(loc = 8) #create legend and label graph axes
    ax.set(title = 'M31 Velocity Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Velocity ($ms^{-1}$))') #label graph axes
    
    """M33 Rotation Curve"""
    #find the velocity profiles of M33 halo and disk
    M33Hv = M33.CircularVelocity(1,r) #vel. profile M33 Halo
    M33Dv = M33.CircularVelocity(2,r) #vel. profile M33 Disk
    M33Tv = M33.CircularVelocityTotal(r) #Total vel. profile
    Hernv = [] #empty array for Herquist vel. profile
    for radius in r: #vel. according to hernquist profile at each radius
        M33Hern = M33.HernquistVCirc(radius,scaleM33*u.kpc,MassM33H)
        Hernv.append(M33Hern*u.s/u.km)
    
    #plot M31 velocity profiles on log graphs
    fig,ax = plt.subplots(figsize = (10,10))
    ax.scatter(r,M33Hv,label = 'Halo Profile')
    ax.scatter(r,M33Dv,label = 'Disk Profile')
    ax.scatter(r,M33Tv,label = 'Total Profile')
    ax.plot(r,Hernv,linestyle = 'dashed',label = 'Hernquist Profile',color = 'red')
    legend=ax.legend() #create legend
    ax.set(title = 'M33 Velocity Profiles', xlabel = 'Radius (kpc)', ylabel = 'Log(Velocity ($ms^{-1}$))') #label graph axes