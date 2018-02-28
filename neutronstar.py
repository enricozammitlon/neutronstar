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
rHat = [[]]
rhoHat = [[]]
mHat = [[]]

def initConstants(rhoS):
    i=0
    rZero=[]
    mZero=[]
    check=[]
    for rhos in rhoS:
        rZero.append((10**-19)*(1/(np.sqrt(rhos*G*4*np.pi)))) #with conversion factor for 10km
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))#with convertion to terms of solar masses.
        check.append((G*4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2)))/(1/(np.sqrt(rhos*G*4*np.pi)))) #check according to chapter should be 1.
        rhoHat[i].append(rhoS[i]/rZero[i])#?
        i = i + 1
    return rZero,mZero,check,rhoHat
print(initConstants(rhoS))


#Im not sure these should return relatible numbers, but the units seem to check out. Mass is in terms of energy and density is energy density.

#solve equations

#function for New pressure to new density

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*n**2.54 +n*MNEUTRON
    return density

def intialPressure(rhoS):
#how do we intialise pressure boio? 
    
# deriv equations for rk4
    
#PLEASE IGNORE/CHANGE BELOW ITS ALL TO GET AN IDEA, its coding crime.
    
def Pderiv(mHat,rHat,densityHat):
    return (-mHat*densityHat/(rHat**2))
#do these need to be functions of can we just put the calcuations into the rk4 arguments? desnityhat is the same as rho hat, that is the same as the dimensionless energry density.
def mderiv(rHat,densityHat):
    return ((rHat**2)*densityHat)
    
#rk4 to give new mass and then pressure.
#rk4 takes y, the diferential and the x value, and aproximates the next point according the the diferential equation, it will spit out y1 and x1 from y0 and x0. 
def rk4(y,dy,x,h):
    k1=dy(y,x)
    k2=dy(y+h/2*k1,x+h/2)
    k3=dy(y+h/2*k2,x+h/2)
    k4=dy(y+h*k3,x+h)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

h = 1

for i in range(0,1000):#because i want to see th result first as we dont know what limit we are looking for 
    (r,m)=rk4(0,mderiv,0,h)#use rk4 for dm/dr, got get m at the next r
    rHat[0[i]].append(r) #append values found for one step, not sure if we need this. but we need to be able to distingush between the diferent rs and ms.
    mHat[0[i]].append(m)
    (r,p)=rk4() # then rk4 for dP/dr m found previous with  to give us pressure at the next r? 
    #do we use the new m with the orginal r? it would not make sense to use the new m and the r it gave us because then we would be taking 2 steps in r for one in pressure.
   #then use pressure to density function to find the new density at the new r.
    #then repeat steps printing P at each step so we can see what happens.







    
