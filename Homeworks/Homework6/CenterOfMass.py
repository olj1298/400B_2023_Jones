#Homework 4
#Center of Mass Position and Velocity
#Solutions: G.Besla, R. Li, H. Foote
#olivia jones 2.9.23
#import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
class CenterOfMass: #Class to define COM position and velocity properties of a given galaxy #and simulation snapshot

    def __init__(self, filename, ptype):
        '''Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            Inputs:
            :filename(string): snapshot file
            :ptype(integer): 1, 2, or 3. particle type to use for COM calculations'''
        
        self.time, self.total, self.data = Read(filename) #read data in the given file using Read                                                                                     
        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)
        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index] #data in m column from file. 
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index] #data in x column from file.  
        self.y = self.data['y'][self.index] #data in y column from file.  
        self.z = self.data['z'][self.index] #data in z column from file.  
        self.vx = self.data['vx'][self.index] #data in vx column from file.  
        self.vy = self.data['vy'][self.index] #data in vy column from file.  
        self.vz = self.data['vz'][self.index] #data in vz column from file.  
    
    def COMdefine(self,a,b,c,m):
        '''Method to compute the COM of a generic vector quantity by direct weighted averaging.
            Inputs:
                :self(function within class):
                :a(float or array of floats): first vector component
                :b(float or array of floats): second vector component
                :c(float or array of floats): third vector component
                :m(float or array of floats): particle masses
            Returns:
                ::a_com(float) first component on the COM vector
                ::b_com(float) second component on the COM vector
                ::c_com(float) third component on the COM vector'''
        #write your own code to compute the generic COM using Eq. 1 in the homework instructions
        #xcomponent Center of mass
        a_com = np.sum(a*m) / np.sum(m)
        #ycomponent Center of mass
        b_com = np.sum(b*m) / np.sum(m)
        #zcomponent Center of mass
        c_com = np.sum(c*m) / np.sum(m)
        
        return a_com,b_com,c_com #return the 3 components separately
    
    def COM_P(self, delta,volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.
        Inputs:
            :self(function within class):
            :delta(float): error tolerance in kpc. Default is 0.1 kpc
            :volDec(float): amount rmax is decreased
        Returns:
            ::p_COM(array of astropy.quantity): 3D position of the center of mass in kpc'''                                                                     
        #Center of Mass Position                                                                                                                                                                      
        #Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        #compute the magnitude of the COM position vector.
        r_COM = np.sqrt(x_COM**2+y_COM**2+z_COM**2)
        
        #iterative process to determine the center of mass                                                            
        #change reference frame to COM frame                                                                          
        #compute the difference between particle coordinates and the first guess at COM position
        x_new = self.x-x_COM #diff in x
        y_new = self.y-y_COM #diff in y
        z_new = self.z-z_COM #diff in z
        r_new = np.sqrt(x_new**2+y_new**2+z_new**2) #magnitude of COM position vector

        #find the max 3D distance of all particles from the guessed COM                                               
        #will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        #pick an initial value for the change in COM position                                                      
        #between the first guess above and the new one computed from half that volume
        #it should be larger than the input tolerance (delta) initially
        change = 1000.0

        #start iterative process to determine center of mass position                                                 
        #delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            #select all particles within the reduced radius (starting from original x,y,z, m)
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2] #x at new radius
            y2 = self.y[index2] #y at new radius
            z2 = self.z[index2] #z at new radius
            m2 = self.m[index2] #mass at new radius

            #Refined COM position:                                                                                    
            #compute the center of mass position using                                                                
            #the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
            #compute the new 3D COM position
            r_COM2 = np.sqrt(x_COM2**2+y_COM**2+z_COM**2) #magnitude of position vector

            #determine the difference between the previous center of mass position and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            #uncomment the following line if you want to check this                                                                                               
            #print("CHANGE = ", change)                                                                                     

            #Before loop continues, reset : r_max, particle separations and COM                                        

            #reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec
            #check this.                                                                                              
            #print("maxR", r_max)                                                                                      

            #Change the frame of reference to the newly computed COM.                                                 
            #subtract the new COM
            x_new = x2 - x_COM2 #diff in x position
            y_new = y2 - y_COM2 #diff in y position
            z_new = z2 - z_COM2 #diff in z position
            r_new = np.sqrt(x_new**2+y_new**2+z_new**2) #magnitude of position vector

            #set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            #create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])
        
        #round to two decimal places
        x_COM = np.round(x_COM,2)    
        y_COM = np.round(y_COM,2) 
        z_COM = np.round(z_COM,2) 

        #set the correct units using astropy and round all values and then return the COM positon vector
        return np.array([x_COM,y_COM,z_COM])*u.kpc 
        
    def COM_V(self, x_COM, y_COM, z_COM):
        '''Method to compute the center of mass velocity based on the center of mass
        position.
        Inputs:
            :self(function within class): 
            :x_COM(astropy quantity): The x component of the center of mass in kpc
            :y_COM(astropy quantity): The y component of the center of mass in kpc
            :z_COM(astropy quantity): The z component of the center of mass in kpc    
        RETURNS:
            ::v_COM(array of astropy quantity): 3D velocity of the center of mass in km/s'''
        #the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc #15kpc mask given in homework

        #determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        #absolute values of position vector from difference
        xV = np.abs(self.x*u.kpc-x_COM) 
        yV = np.abs(self.y*u.kpc-y_COM)
        zV = np.abs(self.z*u.kpc-z_COM)
        rV = np.sqrt(xV**2+yV**2+zV**2) #magnitude of position vector
        
        #determine the index for those particles within the max radius
        indexV = np.where(rV < rv_max) #radius less than max
        
        #determine the velocity and mass of those particles within the mass radius
        vx_new = self.vx[indexV] #velocity in x
        vy_new = self.vy[indexV] #velocity in y
        vz_new = self.vz[indexV] #velocity in z
        m_new =  self.m[indexV] #mass of particles inside radius
        
        #compute the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new,vy_new,vz_new,m_new)
        
        #create an array to store the COM velocity
        v_COM = np.array([vx_COM,vy_COM,vz_COM])

        #return the COM vector set the correct units using astropy round all values                                                                                        
        v_COM = np.round(v_COM,2)*u.km/u.s
        
        return v_COM

if __name__ == '__main__' : 

    #Create a Center of mass object for the MW, M31 and M33
    #below is an example of using the class for MW
    #directory address for Galaxy files
    MWfile = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework4/MW_000.txt"
    M31file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework4/M31_000.txt"
    M33file = "C:/Users/orang/Downloads/400b/400B_2023_Jones/Homeworks/Homework4/M33_000.txt"
    #call and round galaxy mass particle data
    MW_COM = CenterOfMass(MWfile, 2)
    M31_COM = CenterOfMass(M31file, 2)
    M33_COM = CenterOfMass(M33file, 2)
    #below gives you an example of calling the class's functions
    #MW:store the position and velocity COM
    MW_COM_p = MW_COM.COM_P(0.1)
    MW_pmag = np.sqrt(np.sum(np.square(MW_COM_p)))
    print(f"MW COM position components: {MW_COM_p}")
    print(f"MW COM position mag.: {np.round(MW_pmag,2)}")
    MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
    print(f"MW COM vel. components: {MW_COM_v}")
    MW_vmag = np.sqrt(np.sum(np.square(MW_COM_v)))
    print(f"MW COM vel. magnitude: {np.round(MW_vmag,2)}")
    
    #M31:store the position and velocity COM
    M31_COM_p = M31_COM.COM_P(0.1)
    M31_pmag = np.sqrt(np.sum(np.square(M31_COM_p)))
    print(f"M31 COM position components: {M31_COM_p}")
    print(f"M31 COM position mag.:{np.round(M31_pmag,2)}")
    M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
    print(f"M31 COM vel. components: {M31_COM_v}")
    M31_vmag = np.sqrt(np.sum(np.square(M31_COM_v)))
    print(f"M31 COM vel. magnitude: {np.round(M31_vmag,2)}")
    
    #M33:store the position and velocity COM
    M33_COM_p = M33_COM.COM_P(0.1)
    M33_pmag = np.sqrt(np.sum(np.square(M33_COM_p)))
    print(f"M33 COM position components: {M33_COM_p}")
    print(f"M33 COM position mag.: {np.round(M33_pmag,2)}")
    M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
    print(f"M33 COM vel. components: {M33_COM_v}")
    M33_vmag = np.sqrt(np.sum(np.square(M33_COM_v)))
    print(f"M33 COM vel. magnitude: {np.round(M33_vmag,2)}")
    
    #separation between MW and M31
    p_sep = MW_COM_p - M31_COM_p
    v_sep = MW_COM_v - M31_COM_v
    #magnitude of separation
    p_mag = np.sqrt(np.sum(np.square(p_sep)))
    v_mag = np.sqrt(np.sum(np.square(v_sep)))
    print(f"MW and M31 separation: {np.round(p_mag, 3)}")
    print(f"MW and M31 vel. difference: {np.round(v_mag, 3)}")
    #separation between M33 and M31
    p_sep2 = M33_COM_p - M31_COM_p
    v_sep2 = M33_COM_v - M31_COM_v
    #magnitude of separation
    p_mag2 = np.sqrt(np.sum(np.square(p_sep2)))
    v_mag2 = np.sqrt(np.sum(np.square(v_sep2)))
    print(f"M33 and M31 separation: {np.round(p_mag2, 3)}")
    print(f"M33 and M31 vel. difference: {np.round(v_mag2, 3)}")

    #Q4
    print(f"The iterative process was important because the position of the galaxies are moving closer together.")
    print(f"As the galaxies move closer together, the center of mass will change over time as well.")