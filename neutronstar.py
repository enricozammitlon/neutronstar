import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

h=10.0

c=299792458 # m/s
G=(197.327*6.67259*(10**-45)) # not checked
#G=6.67e-11
MSOLAR = 1.1157467*(10**60)/(pow(c,2))#Mev - Energy
MNEUTRON = 938.926/(pow(c,2)) # Mev
HBARC = 197.327 # MevFm
#HBARC=6.58e-7
rhoS=np.arange(19600,25000,1) #MeVFm-3 - Energy desnity from Kgm-3
#do we want rhos in kg and then another array where rhos are in energy density using this converstion?
rHat = [[h] for i in range(len(rhoS)) ]
rhoHat = [[1] for i in range(len(rhoS))]
mHat = [[0] for i in range(len(rhoS))]
pHat = []
rZero=[]
mZero=[]
"""
def f(n):
    return (rhoS[0] - 236*n**(2.54) - n*MNEUTRON)
root = optimize.newton(f,50)
print(root)

root = optimize.newton(f,500)
print(root/rhoS[0])
"""
def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*(n**2.54) +n*MNEUTRON
    return density

def DensityToPressure(rho):
    n=optimize.newton(lambda x: rho-236*(x**(2.54))-x*MNEUTRON,50)
    P = 363.44 * (n**2.54)
    print ("this is p from d2p %6.20f" %P)
    return P

'''
def PressureToDensity1(P):
    rho = pow((3*P*pow(MNEUTRON,4/3)*pow(3*pow(np.pi,2),-1/3))/(HBARC),3/4)
    return rho

def DensityToPressure1(rho):
    P=(HBARC*pow(MNEUTRON,(-4/3))*pow(3*pow(np.pi,2),(1/3))*pow(rho,(4/3))/3)
    return P
'''
def xderiv(phat,rHat,densityHat,mHat):
    if(rHat==0.):
        return 0
    else:
        return (-mHat*densityHat)/pow(rHat,2)
    #relderive
def Pderiv(pHat,rHat,densityHat,mHat):
    if(rHat==0.):
        return 0
    else:
        return(((pHat+densityHat)*(pow(rHat,3)*pHat+mHat))/(pow(rHat,2)-2*mHat*rHat))
    

def mderiv(mHat,rHat,densityHat,b):
    return (pow(rHat,2)*densityHat)

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
        pHat.append([DensityToPressure(rhos)])
    return rZero,mZero,check

print(initConstants(rhoS))

for j in range(0,len(rhoS)):
    #print("FOR DENSITY CENTRAL %f"%(rhoS[j]))
    for i in range(1,5000000):
        #print("Mass:    \t%2.6e"%(mHat[j][i-1]*mZero[j]))
        #print("Pressure:\t%2.6e"%(pHat[j][i-1]))
        #print("New Density:\t%f"%(rhoHat[j][i-1]))
        #print("New R :\t%f"%(PressureToDensity(pHat[j][i-1])))
        (ri,mi)=rk4(mHat[j][i-1], mderiv , rHat[j][i-1], h, rhoHat[j][i-1],mHat[j][i-1])
        (ri,Pi)=rk4(pHat[j][i-1], Pderiv , rHat[j][i-1], h , rhoHat[j][i-1],mHat[j][i-1])
        if(Pi<=0):
            #print("Clause reached")
            break
        pHat[j].append(Pi)
        rHat[j].append(ri)
        mHat[j].append(mi)
        rhoHat[j].append(PressureToDensity(Pi))


fig = plt.figure()
x=[rHat[i][-1]*rZero[i] for i in range(len(rHat))]
y1=[mHat[i][-1]*mZero[i] for i in range(len(mHat))]
plt.xlabel("radius / 10 km")
plt.ylabel("mass / M0")
plt.scatter(x,y1)


fig2 = plt.figure()
y2=[rhoHat[i][0]*rhoS[i] for i in range(len(pHat))]
plt.xlabel("radius / 10 km")
plt.ylabel("density")
plt.scatter(x,y2)

plt.show()

"""
fig = plt.figure()
x=[i*rZero[0] for i in rHat[0]]
y1=[i*mZero[0] for i in mHat[0]]
plt.xlabel("radius /10km")
plt.ylabel("SOLAR MASSES")

plt.plot(x,y1)
plt.show()
"""
