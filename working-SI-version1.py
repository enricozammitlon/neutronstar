# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 17:04:08 2018

@author: cdsch
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

c=299792458
G = 6.67*10**-11
hbar = 1.055*10**(-34)
MNEUTRON = 1.67*10**(-27)

MSolarFinal =[]
RFinal =[]

family = np.arange(20, (2*10**3),50)
for rhoIn in family:

  h = 10.0
  r =[h]
  m =[0]
  MSOL = [0]
  rho = [(rhoIn*10**15)]
  P = []
    
    
  def rk4(y,dy,x,h,rho,m):
    k1=dy(y,x,rho,m)
    k2=dy(y+h/2*k1,x+h/2,rho,m)
    k3=dy(y+h/2*k2,x+h/2,rho,m)
    k4=dy(y+h*k3,x+h,rho,m)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

  def DensityToPressure1(rho):
    P = ((pow(hbar,2)*pow(3*pow(np.pi,2),2/3)*pow(rho,5/3))/(5*pow(MNEUTRON,8/3)))
    return P

  def PressureToDensity1(P):
    rho = pow(((P*5*(pow(MNEUTRON,8/3)))/(pow(hbar,2)*(pow(3*(pow(np.pi,2)),2/3)))),3/5)
    return rho

  def PressureToDensity2(P):
    rho = pow((3*P*pow(MNEUTRON,4/3)*pow(3*pow(np.pi,2),-1/3))/(hbar*c),3/4)
    return rho
    
  def DensityToPressure2(rho):
    P=(hbar*c*pow(MNEUTRON,(-4/3))*pow(3*pow(np.pi,2),(1/3))*pow(rho,(4/3))/3)
    return P

  def PressureToDensity3(Pressure):
    MNEUTRON2=989
    p=(Pressure*10**-45)/(1.602e-13)
    #P = 363.44 * n**2.54
    n = (p/363.44)**(1/2.54)
    density = 236*(n**2.54) +n*MNEUTRON2
    rho=(density*1.602e-13)/(1e-45*c**2)
    #print ("this is rho from p2d %2.6e" %rho)
    return rho

  def DensityToPressure3(rho):
    MNEUTRON2=989
    rho2=(rho*10**-45*c**2)/(1.602*10**-13)
    print("rho2 is %2.2f"%rho2)
    n=optimize.newton(lambda x: rho2-236*(x**(2.54))-x*MNEUTRON2,0.5)
    P = 363.44 * (n**2.54)
    P2 = (P*1.602*10**-13)/(10**-45)
    print ("this is p from d2p %2.6e" %P2)
    return P2

  def Xderiv(p,r,density,m):
    if(r==0):
        return 0
    else:
        return ((-G*m*density)/pow(r,2))

  def Pderiv(P,r,rho,m): # working relivitistic 
    if m == 0:
        return 0
    else:
        a = (1+(P/(rho*pow(c,2))))*(1+(4*np.pi*pow(r,3)*P)/(m*pow(c,2)))*pow((1-(2*G*m/(r*pow(c,2)))),-1)
        b = (-G*m*rho/pow(r,2))
        d=a*b
        return d
  
  def Mderiv(M,r,rho,b):
      return (4*rho*np.pi*pow(r,2))
  
  P.append(DensityToPressure3(rho[0]))
  print("IDEALLY2: %2.6e"%(DensityToPressure2(rho[0])))
  #P.append(500)
  for i in range(1,100000):
    """
    #if rho[i-1] >  3.8*10*17:
    if rho[i-1] > 2e17:
      print("Mass:    \t%2.6e"%(m[i-1]))
      print("Pressure:\t%2.6e"%(P[i-1]))
      print("New Density:\t%e"%(rho[i-1]))
      print("New R :\t%f"%(r[i-1]))
      (r1,m1)= rk4(m[i-1], Mderiv ,r[i-1],h,rho[i-1],1)
      (r1,P1) = rk4(P[i-1], Pderiv ,r[i-1],h,rho[i-1],m[i-1])
      if(P1<=0):
          MSolarFinal.append((m[i-1]/(1.989*10**30)))
          RFinal.append(r[i-1]/1000)
          break
      P.append(P1)
      r.append(r1)
      m.append(m1)
      MSOL.append(m1/(1.989*10**30))
      rho.append(PressureToDensity2(P1))
    else: 
      print("Mass:    \t%2.6e"%(m[i-1]))
      print("Pressure:\t%2.6e"%(P[i-1]))
      print("New Density:\t%e"%(rho[i-1]))
      print("New R :\t%f"%(r[i-1]))
      (r1,m1)= rk4(m[i-1], Mderiv ,r[i-1],h,rho[i-1],1)
      (r1,P1) = rk4(P[i-1], Pderiv ,r[i-1],h,rho[i-1],m[i-1])
      if(P1<=0):
          MSolarFinal.append((m[i-1]/(1.989*10**30)))
          RFinal.append(r[i-1]/1000)
          break
      P.append(P1)
      r.append(r1)
      m.append(m1)
      MSOL.append(m1/(1.989*10**30))
      rho.append(PressureToDensity1(P1))
    """
    #print("Mass:    \t%2.6e"%(m[i-1]))
    #print("Pressure:\t%2.6e"%(P[i-1]))
    #print("New Density:\t%e"%(rho[i-1]))
    #print("New R :\t%f"%(r[i-1]))
    (r1,m1)= rk4(m[i-1], Mderiv ,r[i-1],h,rho[i-1],1)
    (r1,P1) = rk4(P[i-1], Pderiv ,r[i-1],h,rho[i-1],m[i-1])
    if(P1<=0):
      MSolarFinal.append((m[i-1]/(1.989*10**30)))
      RFinal.append(r[i-1]/1000)
      break
    P.append(P1)
    r.append(r1)
    m.append(m1)
    MSOL.append(m1/(1.989*10**30))
    #print("IDEALLY: %2.6e"%(PressureToDensity1(P1)))
    rho.append(PressureToDensity3(P1))

print (MSolarFinal, RFinal)

fig = plt.figure()
x1=RFinal
y=MSolarFinal
y2=RFinal
x2=(family*10**15)
plt.xlabel("radius KM")
plt.ylabel("SOLAR MASSES")
plt.plot(x1, y, "o")


fig2=plt.figure()
plt.plot(x2, y, "o")
plt.xlabel("central density kg/m^3")
plt.ylabel("SOLAR MASSES")
plt.show()

fig3=plt.figure()
plt.plot(x2, y2, "o")
plt.xlabel("central desnisty kg/m^3")
plt.ylabel("R")
plt.show()

    
