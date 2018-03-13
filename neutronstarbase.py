import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


G= 1.327400357356e-42 # not checked
den = 5.586592179e-16 #converts kg - mev/c^2fm^3
MSOLAR = 1.989e30/(1.79e-30)#Mev/c^2
MNEUTRON = 938.926 # Mev/c^2
HBAR = 4.135667662e-22# MevS


rhoS=(5e17)*den #MeVFm-3c^-2 - Energy desnity from Kgm-3
R0 = pow(4*np.pi*G*rhoS,-1/2)*1e-19
M0 = (4*np.pi*rhoS)/(pow(4*np.pi*rhoS*G,3/2)*MSOLAR)
alpha=3.485e-4 * pow(rhoS,2/3) # constant of proportionality for EoS - non rel
beta = 1 # constant of proportionality for the rel.

print(rhoS)
#do we want rhos in kg and then another array where rhos are in energy density using this converstion?
h=0.001
rHat = [h]
rhoHat = [1]
mHat = [0]
pHat = []
rZero=[h*R0]
mZero=[0]
pZero=[]

def f(n):
    return (rhoS - 236*n**(2.54) - n*MNEUTRON)
root = optimize.newton(f,500)
print(root)
'''
Returns P hat
input is rho
'''
'''
def DensityToPressure(rhoHat):
        P = alpha*pow(rhoHat,5/3)
        return P
'''   
'''
Return r hat
input is P hat
'''
'''
def PressureToDensity(P):
    return pow(P/alpha,3/5)
'''

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*(n**2.54) +n*MNEUTRON
    return density

def DensityToPressure(n):
    P = 363.44 * n**2.54
    print ("this is p from d2p %6.20f" %n)
    return P


def Pderiv(phat,rHat,densityHat,mHat):
    '''
    if(rHat==0.):
        return 0
    else:
        '''
    return (-mHat*densityHat/pow(rHat,2))

def PderivRel(pHat,rHat,denistyHat,mHat):
    return(((pHat+rhoHat)*(pow(rHat,3))/(pow(rHat - 2*mHat*rHat))))

def mderiv(mHat,rHat,densityHat,b):
    return pow(rHat,2)*densityHat

pHat.append(DensityToPressure(rhoS)) #intialise
pZero.append(pHat[0]*1.6e32*rhoS) #intialise, 1.6e32 changes mev to pascal

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
'''
def initConstants(rhoS):
    check=[]
    for rhos in rhoS:
        rZero.append((10**-19)*(1/(np.sqrt(rhos*G*4*np.pi)))) #with conversion factor for 10km
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))#with convertion to terms of solar masses.
        check.append((G*4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2)))/(1/(np.sqrt(rhos*G*4*np.pi)))) #check according to chapter should be 1.
        pHat.append([DensityToPressure(root/rhoS[0])])
    return rZero,mZero,check
'''

for i in range(1,100000000):
    
   # if:
        
    '''print("Mass:    \t%2.6e"%(mHat[j][i-1]*mZero[j]))
    print("Pressure:\t%2.6e"%(pHat[j][i-1]))
    print("New R:\t%f"%(rHat[j][i-1]*rZero[j]))
    print("New Density :\t%f"%(PressureToDensity(pHat[j][i-1])))'''
    (ri,mi)=rk4(mHat[i-1], mderiv , rHat[i-1], h, rhoHat[i-1],mHat[i-1])
    rHat.append(ri)
    mHat.append(mi)
    rZero.append(ri*R0)
    mZero.append(mi*M0)
    (ri,Pi)=rk4(pHat[i-1], Pderiv , rHat[i-1], h , rhoHat[i-1],mHat[i-1])
    pHat.append(Pi)
    pZero.append(Pi*1.6e32*rhoS)
    
    if(Pi<=0):
        break
    
    rhoHat.append(PressureToDensity(Pi))
      
        

fig = plt.figure()
x= rZero
y1=mZero
y2=pZero
plt.xlabel("radius 10Km")
plt.ylabel("mass Solar")

plt.plot(x,y1)

fig2=plt.figure()
plt.plot(x,y2)
plt.xlabel("radius 10Km")
plt.ylabel("pressure Pa")
plt.show()
