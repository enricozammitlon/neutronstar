
import numpy as np

MNEUTRON = 938.926 # Mev
rho = [2]
P=[10000]
m=[4]
r=[1]

def rk4(y,dy,x,h,rhos):
    k1=dy(y,x,rhos)
    k2=dy(y+h/2*k1,x+h/2,rhos)
    k3=dy(y+h/2*k2,x+h/2,rhos)
    k4=dy(y+h*k3,x+h,rhos)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*n**2.54 +n*MNEUTRON
    return density

def Pderiv(mHat,rHat,densityHat):
    return (-mHat*densityHat/(rHat**2))

def mderiv(m,rHat,densityHat):
    return ((rHat**2)*densityHat)



h=0.1

for i in range(1,50):
    print("Mass%6.2e : "%(m[i-1]))
    (ri,mi)=rk4(m[i-1], mderiv , r[i-1], h, rho[i-1])
    r.append(ri)
    m.append(mi)
    (ri,Pi)=rk4(P[i-1], Pderiv , r[i-1],h,rho[i-1])
    P.append(Pi)
    rho.append(PressureToDensity(Pi))
    print("Pressure%6.2e :"%Pi)



    #was working on this as practisse
