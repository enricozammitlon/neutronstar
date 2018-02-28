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

# deriv equations for rk4
    
#PLEASE IGNORE BELOW ITS ALL TO GET AN IDEA, its coding crime.
    
def Pderiv(mHat,rHat,densityHat):
    return (-mHat*densityHat/(rHat**2))

def mderiv(rHat,densityHat):
    return ((rHat**2)*densityHat)
    
#rk4 to give new mass and then pressure. ( each new pressure )
    
def rk4(y,dy,x,h):
    k1=dy(y,x)
    k2=dy(y+h/2*k1,x+h/2)
    k3=dy(y+h/2*k2,x+h/2)
    k4=dy(y+h*k3,x+h)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

h = 1

for i in range(0,1000):#because i want to see th result firstas we dont know what limit we are looking for 
    (r,m)=rk4(0,mderiv,0,h)#ive gone blank
    rHat[0[i]].append(r)
    mHat[0[i]].append(m)
    (p,r)=rk4()
    #i wnt to use the resultant m here with the original r? in the nex 






    
