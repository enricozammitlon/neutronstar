# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 17:04:08 2018

@author: cdsch
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

#Important constants
c=299792458 #Speed of light
G = 6.67*10**-11  #Gravitational Constant
hbar = 1.055*10**(-34) #HBAR
MNEUTRON = 1.67*10**(-27) #Mass of a neutron

"""
Will hold the final mass of the star in solar masses.
Each sub array inside the 2D array holds a different method, namely,the
D2P1 initialisation with P2D2 and P2D3 in the first one, and D2P3 with P2D3 in
the second one. The same applies for the RFinal 2D array but holding radius.
"""
MSolarFinal =[[],[]]
RFinal =[[],[]]
#This is used to initialise the boundary for central densities
family = np.arange(1.69, (1.09*10**2),0.5) #To be multipled x10^15
#This dictionary serves to then label the graphs appropriately, make use of it!
methodNames={0:'Other',1:'Bethe And Johnson'}
#Our RK4 integration numerical method method
def rk4(y,dy,x,h,rho,m):
    k1=dy(y,x,rho,m)
    k2=dy(y+h/2*k1,x+h/2,rho,m)
    k3=dy(y+h/2*k2,x+h/2,rho,m)
    k4=dy(y+h*k3,x+h,rho,m)
    y=y+h*(k1+2*k2+2*k3+k4)/6
    x=x+h
    return (x,y)

#To be used when?
def DensityToPressure1(rho):
  P = ((pow(hbar,2)*pow(3*pow(np.pi,2),2/3)*pow(rho,5/3))/(5*pow(MNEUTRON,8/3)))
  return P

#Used in lower density ranges, with D2P2
def PressureToDensity1(P):
  rho = pow(((P*5*(pow(MNEUTRON,8/3)))/(pow(hbar,2)*(pow(3*(pow(np.pi,2)),2/3)))),3/5)
  return rho

#Used in higher density ranges, with D2P2
def PressureToDensity2(P):
  rho = pow((3*P*pow(MNEUTRON,4/3)*pow(3*pow(np.pi,2),-1/3))/(hbar*c),3/4)
  return rho

#To be used with P2D1 and P2D2
def DensityToPressure2(rho):
  P=(hbar*c*pow(MNEUTRON,(-4/3))*pow(3*pow(np.pi,2),(1/3))*pow(rho,(4/3))/3)
  return P

#To be used with D2P3
def PressureToDensity3(Pressure):
  MNEUTRON2=989
  p=(Pressure*10**-45)/(1.602e-13)
  #P = 363.44 * n**2.54
  n = (p/363.44)**(1/2.54)
  density = 236*(n**2.54) +n*MNEUTRON2
  rho=(density*1.602e-13)/(1e-45*c**2)
  return rho

#To be used with P2D3 and is EOS from the chapter given.
def DensityToPressure3(rho):
  MNEUTRON2=989
  rho2=(rho*10**-45*c**2)/(1.602*10**-13)
  n=optimize.newton(lambda x: rho2-236*(x**(2.54))-x*MNEUTRON2,0.5)
  P = 363.44 * (n**2.54)
  P2 = (P*1.602*10**-13)/(10**-45)
  return P2

#This is the classical derivative for pressure.
#!!!!Rename to Pderiv if you want to use this instead of the TOV!!!!
def Xderiv(p,r,density,m):
  if(r==0):
    return 0
  else:
    return ((-G*m*density)/pow(r,2))

#This is the relativistic form of the derivative aka TOV
def Pderiv(P,r,rho,m):
  if m == 0:
    return 0
  else:
    a = (1+(P/(rho*pow(c,2))))*(1+(4*np.pi*pow(r,3)*P)/(m*pow(c,2)))*pow((1-(2*G*m/(r*pow(c,2)))),-1)
    b = (-G*m*rho/pow(r,2))
    d=a*b
    return d

#The mass derivative
def Mderiv(M,r,rho,b):
  return (4*rho*np.pi*pow(r,2))

#Extrapolation technique to be used towards the last 2 points before pressure
#Becomes 0
def Extrap(p1,p2,m1,m2,r1,r2):
  grad1 = (p1-p2)/(r1-r2) #Take gradient of pressure
  c = P1 - grad1*r1 #This is the linear form y=mx+c
  #Use newton-raphson numerical analysis to get the radius at P=0
  RadExtrap =optimize.newton(lambda R: +grad1*R +c, r1)
  grad2 = (m1-m2)/(r1-r2) #Mass gradient
  d = m1 - grad2*r1 #Calculate offset of mass using y=mx+c again
  MassExtrap = grad2*RadExtrap + d #Extrapolate to where P=0 ( R(P=0))
  return RadExtrap,MassExtrap

#Number of combinations of P2D and D2P used. Recall (n+1) since 0 is included
methods=2

for currentmethod in range(methods): #For each method combination
    for rhoIn in family: #For each central density

      h = 10.0 #Step size in meters
      r =[h] #Will hold the radius as it grows
      m =[0] #Will hold the mass as it grows
      rho = [(rhoIn*10**17)] #Will hold the density as it decreases
      P = [] #Will hold the pressure as it decreases
      breakpoint=1e17 #For method 0 a breakpoint between high and low densities is needed

      #Initialise the pressure using the appropriate D2P(current_central_density)
      if currentmethod==0 :
          P.append(DensityToPressure2(rho[0]))
      if currentmethod==1:
          P.append(DensityToPressure3(rho[0]))

      for i in range(1,100000): #Arbitrary number of iterations to ensure enough
          #print("Mass:    \t%2.6e"%(m[i-1]))
          #print("Pressure:\t%2.6e"%(P[i-1]))
          #print("New Density:\t%e"%(rho[i-1]))
          #print("New R :\t%f"%(r[i-1]))

          #Use the previous mass,radius and density to calculate the next one
          (r1,m1)= rk4(m[i-1], Mderiv ,r[i-1],h,rho[i-1],1)
          #Use the previous mass,radius,density and pressure to calculate the next one
          (r1,P1) = rk4(P[i-1], Pderiv ,r[i-1],h,rho[i-1],m[i-1])
          if(P1<=0):#Boundary condition to know when to stop since edge is reached
              #Extrapolate the pressure and mass to where P=0
              (RF,MF) = Extrap(P[i-1],P[i-2],m[i-1],m[i-2],r[i-1],r[i-2])
              print("Final Mass: %f"%MF)
              print("Final Radius: %f"%RF)
              #Append the final star mass and radius for the given method
              MSolarFinal[currentmethod].append((MF/(1.989*10**30)))
              RFinal[currentmethod].append(RF/1000)
              break
          #Append these to be used for the (i+1)th pressure and mass finding
          P.append(P1)
          r.append(r1)
          m.append(m1)
          #Decide which P2D method to use depending on the current method being investigated
          if(rho[i-1]>breakpoint and currentmethod==0):
              rho.append(PressureToDensity2(P1))
          elif(rho[i-1]<=breakpoint and currentmethod==0):
              rho.append(PressureToDensity1(P1))
          elif(currentmethod==1):
              rho.append(PressureToDensity3(P1))

#A graph of Mass vs Radius for the whole family,for all methods
fig = plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    x2=(family*10**15)
    plt.xlabel("Radius/km")
    plt.ylabel("Solar Masses")
    plt.plot(x1, y, "o",label=methodNames[i])
plt.legend()

#A graph of Mass vs Central Density for the whole family,for all methods
fig2=plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    x2=(family*10**15)
    plt.plot(x2, y, "o",label=methodNames[i])
    plt.xlabel("Central Density/ kgm^-3")
    plt.ylabel("Solar Masses")
plt.legend()
#A graph of Radius vs Central Density for the whole family,for all methods
fig3=plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    x2=(family*10**15)
    plt.plot(x2, y2, "o",label=methodNames[i])
    plt.xlabel("Central Density/ kgm^-3")
    plt.ylabel("Radius/km")
plt.legend()
plt.show()
