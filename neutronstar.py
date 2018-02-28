import sympy as sp
import numpy as np

global G
global MSOLAR
global c
c=299792458 # m/s
G=197.327*6.67259*(10**-45) # not checked
MSOLAR=1.1157467*(10**60) #Mev - Energy
MNEUTRON = 938.926 # Mev
rhoS=[((3.7*10**17)*10**(-45)*(3*10**8)**2)/(1.6*10**(-13))] #MeVFm-3 - Energy desnity from Kgm-3
#do we want rhos in kg and then another array where rhos are in energy density using this converstion? 

def initConstants(rhoS):
    rZero=[]
    mZero=[]
    check=[]
    for rhos in rhoS:
        rZero.append((10**-19)*(1/(np.sqrt(rhos*G*4*np.pi)))) #with conversion factor for 10km
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))#with convertion to terms of solar masses.
        check.append((G*4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2)))/(1/(np.sqrt(rhos*G*4*np.pi)))) #check according to chapter should be 1.
    return rZero,mZero,check
print(initConstants(rhoS))
#Im not sure these should return relatible numbers, but the units seem to check out. Mass is in terms of energy and density is energy density.

#solve equations

#function for New pressure to new density

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*n**2.54 +n*MNEUTRON
    return density

#rk4 to give new pressure. ( each new pressure )
    
