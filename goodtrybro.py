
import numpy as np

MNEUTRON = 938.926 # Mev
rho = [2]
P=[10000]
m=[4]
r=[1]

def rk4m(y,dy,x,h,rhos):
    k1=dy(y,x,rhos)
    k2=dy(y+h/2*k1,x+h/2,rhos)
    k3=dy(y+h/2*k2,x+h/2,rhos)
    k4=dy(y+h*k3,x+h,rhos)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

def rk4P(y,dy,x,h,rhos):
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



run= np.arange(1,10000)
h=0.1

for i in run:
    print(m[i-1])
    (ri,mi)=rk4m(m[i-1], mderiv(m[i-1],r[i-1],rho[i-1]), r[i-1], h, rho[i-1])
    r.append(ri)
    m.append(mi)
    (ri,Pi)=rk4(P[i-1], Pderiv(m[i-1],r[i-1],rho[i-1]), r[i-1],h)
    P.append(Pi)
    rho.appendPressureToDensity(Pi)
    print(Pi)
    
    
    
    #was working on this as practisse
    
    
    
    
    
    
    
    
    
    
    