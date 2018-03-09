import numpy as np
from scipy.optimize import fsolve

C=299792458 # m/s
G=6.674e-11
MSOLAR=1.989e30
MNEUTRON = 1.67e-27
HBAR = 1.05e-34
rho=[2.3e17]
p=[]
m=[]
r=[]

"""
def DensityToPressure(rho,rhoprev):
    pressure= ( (HBAR**2) * (3*(np.pi**2))**(2/3) * rho**(5/3) ) / (5*(MNEUTRON**(8/3)))
    s=fsolve(lambda x : rho-((5*x*(MNEUTRON**(8/3)))/((HBAR**2)*(3*pow(np.pi,2))**(2/3)))**(3/5),rhoprev)
    return pressure

def PressureToDensity(p):
    density = ((5*p*(MNEUTRON**(8/3)))/((HBAR**2)*(3*pow(np.pi,2))**(2/3)))**(3/5)
    return density
"""
#2.04e17
def PressureToDensity(p,rho):
    s=fsolve(lambda x : p-9.58e-28*x**2.54,rho)
    return s

def DensityToPressure(rho):
    p=(9.58e-28)*rho**2.54
    return p

print(PressureToDensity(DensityToPressure(rho[0]),8e-17))
"""
dp = -G*m*rho/c^2 r^2
dr
"""
def pDeriv(p,r,rho,m):
    if(r==0):
        return 0
    return (-G*m*rho)/(r**2 * C**2)

"""
dm  = 4πr^2 * ρ(r) / c^2 .
dr
"""
def mDeriv(m,r,rho,extra):
    return (4*np.pi*(r**2)*rho)/C**2

def rk4(y,dy,x,h,rho,m):
    k1=dy(y,x,rho,m)
    k2=dy(y+h/2*k1,x+h/2,rho,m)
    k3=dy(y+h/2*k2,x+h/2,rho,m)
    k4=dy(y+h*k3,x+h,rho,m)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)
print(DensityToPressure(rho[0]))
p.append(1.6e34)
m.append(0)
r.append(0)

h=0.2

for i in range(1,5):
    print("Mass:    \t%2.6e"%(m[i-1]))
    print("Pressure:\t%2.6e"%(p[i-1]))
    print("Radius  :\t%2.6e"%(r[i-1]))
    print("Density:\t%2.6e"%(rho[i-1]))
    print("_______________________________")
    (ri,mi)=rk4(m[i-1], mDeriv , r[i-1], h, rho[i-1],m[i-1])
    r.append(ri)
    m.append(mi)
    (ri,Pi)=rk4(p[i-1], pDeriv , r[i-1], h , rho[i-1],m[i-1])
    p.append(Pi)
    if(Pi<=0):
        print("CLAUSE REACHED")
        break
    rho.append(PressureToDensity(Pi,rho[i-1]))
