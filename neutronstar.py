import numpy as np

c=299792458 # m/s
G=197.327*6.67259*(10**-45) # not checked
MSOLAR=1.1157467*(10**60) #Mev - Energy
MNEUTRON = 938.926 # Mev
HBARC = 197.327 # MevFm

rhoS=[((3.7*10**17)*10**(-45)*(3*10**8)**2)/(1.6*10**(-13))] #MeVFm-3 - Energy desnity from Kgm-3
#do we want rhos in kg and then another array where rhos are in energy density using this converstion?
rHat = [[0]]
rhoHat = [[1]]
mHat = [[0]]
pHat = []
m=[[0]]
r=[[0]]
p=[]
rZero=[]
mZero=[]
rho=[]

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*(n**2.54) +n*MNEUTRON
    return density

def CentralPressure(rho):#this matches the units required
    PressureCentral= ((HBARC*(pow(3*np.pi,2)*pow(rho,(5/3))))/(5*MNEUTRON)**(8/3))
    return PressureCentral


# deriv equations for rk4

def Pderiv(phat,rhat,densityHat,mhat):
    if(rhat==0.):
        return 0
    else:
        return (-mhat*densityHat/(rhat**2))

def mderiv(mHat,rHat,densityHat,b):
    return ((rHat**2)*densityHat)

#rk4 to give new mass and then pressure.
#rk4 takes y, the diferential and the x value, and aproximates the next point according the the diferential equation, it will spit out y1 and x1 from y0 and x0.
def rk4(y,dy,x,h,rhos,mHat):
    k1=dy(y,x,rhos,mHat)
    k2=dy(y+h/2*k1,x+h/2,rhos,mHat)
    k3=dy(y+h/2*k2,x+h/2,rhos,mHat)
    k4=dy(y+h*k3,x+h,rhos,mHat)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

def initConstants(rhoS):
    check=[]
    for rhos in rhoS:
        rZero.append((10**-19)*(1/(np.sqrt(rhos*G*4*np.pi)))) #with conversion factor for 10km
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))#with convertion to terms of solar masses.
        check.append((G*4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2)))/(1/(np.sqrt(rhos*G*4*np.pi)))) #check according to chapter should be 1.
        pHat.append([CentralPressure(rhos)/rhos])
    return rZero,mZero,check

print(initConstants(rhoS))

h=0.000001
for j in range(0,len(rhoS)):
    print("FOR DENSITY CENTRAL %f"%(rhoHat[j][0]))
    for i in range(1,5):
        print("Mass:    \t%2.6e"%(mHat[j][i-1]/mZero[j]))
        print("Pressure:\t%2.6e"%(pHat[j][i-1]))
        print("New Density:\t%f"%(PressureToDensity(pHat[j][i-1])))
        (ri,mi)=rk4(mHat[j][i-1], mderiv , rHat[j][i-1], h, rhoHat[j][i-1],mHat[j][i-1])
        rHat[j].append(ri)
        mHat[j].append(mi)
        (ri,Pi)=rk4(pHat[j][i-1], Pderiv , rHat[j][i-1], h , rhoHat[j][i-1],mHat[j][i-1])
        pHat[j].append(Pi)
        if(Pi<=0):
            break
        rhoHat[j].append(PressureToDensity(Pi))
