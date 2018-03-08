import numpy as np

from scipy import optimize
import matplotlib.pyplot as plt




c=299792458 # m/s
G=1.67*1088-11 # not checked
MSOLAR=1.989*30 #Mev - Energy
MNEUTRON = 1.67*10**-27 # Mev
#HBARC = 197.327 # MevFm

rhoS=[(3.7*10**15)] #MeVFm-3 - Energy desnity from Kgm-3
#rhoS=[10000]
#do we want rhos in kg and then another array where rhos are in energy density using this converstion?
r = [[0]]
rho = [[rhoS[0]]]
m = [[0]]
p =[]

def f(n):
    return (rhoS[0] - 236*n**(2.54) - n*MNEUTRON)

def IntialPressure(rhoS):
    root = optimize.newton(f,500000000)
    print(" n = %6.10f" %root)
    n = root
    P0 = 363.44 * n**2.54
    return P0

def PressureToDensity(Pressure):
    #P = 363.44 * n**2.54
    n = (Pressure/363.44)**(1/2.54)
    density = 236*(n**2.54) +n*MNEUTRON
    return density

def DensityToPressure(n):
    P = 363.44 * n**2.54
    return P

def Pderiv(p,r,density,m):
    if(r==0.):
        return 0
    else:
        return (-G*m*density/(r**2))
    

def mderiv(m,r,density,b):
    return ((4*np.pi*r**2)*density)

#rk4 to give new mass and then pressure.
#rk4 takes y, the diferential and the x value, and aproximates the next point according the the diferential equation, it will spit out y1 and x1 from y0 and x0.
def rk4(y,dy,x,h,rho,m):
    k1=dy(y,x,rhos,m)
    k2=dy(y+h/2*k1,x+h/2,rhos,m)
    k3=dy(y+h/2*k2,x+h/2,rhos,m)
    k4=dy(y+h*k3,x+h,rho,m)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

for rhos in rhoS:
   p.append([IntialPressure(rho[0])])
    
h=0.000000000001

for j in range(0,len(rhoS)):
    print("FOR DENSITY CENTRAL %f"%(rho[j][0]))
    for i in range(1,5):
        print("Mass:    \t%2.6e"%(m[j][i-1]))
        print("Pressure:\t%2.6e"%(p[j][i-1]))
        print("New Density:\t%f"%(PressureToDensity(p[j][i-1])))
        (ri,mi)=rk4(m[j][i-1], mderiv , r[j][i-1], h, rho[j][i-1],m[j][i-1])
        r[j].append(ri)
        m[j].append(mi)
        (ri,Pi)=rk4(p[j][i-1], Pderiv , r[j][i-1], h , rho[j][i-1],m[j][i-1])
        p[j].append(Pi)
        if(Pi<=0):
            break
        rho[j].append(PressureToDensity(Pi))
       # print(mHat[j][i-1]*mZero[j])
       



fig = plt.figure()
x=r[0]
y1=m[0]
y2=p[0]
plt.xlabel("radius meters")
plt.ylabel("mass kg")

plt.plot(x,y1)

fig2=plt.figure()
plt.plot(x,y2)
plt.xlabel("radius meters")
plt.ylabel("pressure Nm^-2")
plt.show()

