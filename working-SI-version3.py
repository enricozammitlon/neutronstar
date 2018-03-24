# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 17:04:08 2018

@authors: cdsch,enrico
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from matplotlib.patches import Wedge
import os
import math
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
RFError =[[],[]]
MFE =[[],[]]
RFE =[[],[]]
#This is used to initialise the boundary for central densities
family = np.arange(1.69, (1.09*10**2),1) #To be multipled x10^17 using BJ resitrcitions
#This dictionary serves to then label the graphs appropriately, make use of it!
methodNames={0:'Ideal Neutron Degenerate Gas',1:'Bethe And Johnson'}
#Our RK4 integration numerical method method
def rk4(y,dy,x,h,rho,m):
    results=[]
    for i in range(len(y)):#This loops for all the values of h,2h...etc
        k1=dy(y[i],x[i],rho[i],m[i])
        k2=dy(y[i]+h[i]/2*k1,x[i]+h[i]/2,rho[i],m[i])
        k3=dy(y[i]+h[i]/2*k2,x[i]+h[i]/2,rho[i],m[i])
        k4=dy(y[i]+h[i]*k3,x[i]+h[i],rho[i],m[i])
        yf=y[i]+h[i]*(k1+2*k2+2*k3+k4)/6
        xf=x[i]+h[i]
        results.append([xf,yf])
    return results

#our RK5 method for giving a comparitive againt the rk4's accuracy
#( the rk5 is theoretcally more accurate but, numerical errors are larger
#this will be used in disucssion as to how good the rk4 is.
def rk5(y,dy,x,h,rho,m):
    results=[]
    for i in range(len(y)):#This loops for all the values of h,2h...etc
        k1=dy(y[i],x[i],rho[i],m[i])
        k2=h[i]*dy(y[i]+k1/2,x[i]+h[i]/2,rho[i],m[i])
        k3=h[i]*dy(y[i]+(3*k1+k2)/16,x[i]+h[i]/4,rho[i],m[i])
        k4=h[i]*dy(y[i]+k3/2,x[i]+h[i]/2,rho[i],m[i])
        k5=h[i]*dy(y[i]+(-3*k2+6*k3+9*k4)/16, x[i] + 3*h[i]/4,rho[i],m[i])
        k6=h[i]*dy(y[i]+(k1+4*k2+6*k3-12*k4+8*k5)/7, x[i] + h[i],rho[i],m[i])
        yf=y[i]+(7*k1+32*k3+12*k4+32*k5+7*k6)/90
        xf=x[i]+h[i]
        results.append([xf,yf])
    return results

def Euler(y,dy,x,h,rho,m):
    results=[]
    for i in range(len(y)):#This loops for all the values of h,2h...etc
        yf = y[i] + h[i]*dy(y[i],x[i],rho[i],m[i])
        xf = x[i]+h[i]
        results.append([xf,yf])
    return results


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
omega = 2*np.pi*0
#omega is used to investigate the addition of (special [not accelerating]) relavtatistic fictious forces
#set omega to zero to ingore this for the main data set
def Xderiv(p,r,density,m):
  if(r==0):
    return 0
  else:
    u = (-G*m*density)/pow(r,2)
    v = r * omega
    gamma = 1/pow((1-(pow(v,2)/pow(c,2))),1/2)
    s= gamma*density*pow(v,2)/r
    return (u+s)

#This is the relativistic form of the derivative aka TOV


def Pderiv(P,r,rho,m):
  if m == 0:
    return 0
  else:
    a = (1+(P/(rho*pow(c,2))))*(1+(4*np.pi*pow(r,3)*P)/(m*pow(c,2)))*pow((1-(2*G*m/(r*pow(c,2)))),-1)
    b = (-G*m*rho/pow(r,2))
    v = r * omega
    gamma = 1/pow((1-(pow(v,2)/pow(c,2))),1/2)
    s= gamma*rho*pow(v,2)/r
    d=(a*b)+s
    return d

#The mass derivative
def Mderiv(M,r,rho,b):
  return (4*rho*np.pi*pow(r,2))

#Extrapolation technique to be used towards the last 2 points before pressure
#Becomes 0
def Extrap(p1,p2,m1,m2,r1,r2,h):
  grad1 = (p1-p2)/(r1-r2) #Take gradient of pressure
  c = p1 - grad1*r1 #This is the linear form y=mx+c
  #Use newton-raphson numerical analysis to get the radius at P=0
  RadExtrap =optimize.newton(lambda R: +grad1*R +c, r1)
  grad2 = (m1-m2)/(r1-r2) #Mass gradient
  d = m1 - grad2*r1 #Calculate offset of mass using y=mx+c again
  MassExtrap = grad2*RadExtrap + d #Extrapolate to where P=0 ( R(P=0))
  #error in extrapolations for linear is (x1-x2)^2! This will be combined with the rk4 trunkation error prportional to h^5
  RadError = pow((r1/1000-r2/1000),2)
  return RadExtrap,MassExtrap,RadError

#Number of combinations of P2D and D2P used. Recall (n+1) since 0 is included
methods=2
flags=[[],[]]

for currentmethod in range(methods):#For each method combination
    counter=0
    for rhoIn in family: #For each central density
      counter+=1
      os.system('clear')
      print("Method: %d"%currentmethod)
      for x in range(counter):
          if(x%10==0):
              print("-",end='', flush=True)
      print("> %d%% Done"%(counter/109 *100),end='', flush=True)

      h = [10,20] #Step size in meters
      r =[[h[0]],[h[1]]] #Will hold the radius as it grows
      m =[[0],[0]] #Will hold the mass as it grows
      rho = [[(rhoIn*1e17)],[(rhoIn*1e17)]] #Will hold the density as it decreases
      P = [] #Will hold the pressure as it decreases
      breakpoint=2.3e17 #For method 0 a breakpoint between high and low densities is needed
      mrk5 =[0]#further arrays for the rk5 comparitive method, if used.
      Prk5 = []
      boundOne=False
      boundTwo=False
      boundThree=False
      continueSecondStar=True
      #Initialise the pressure using the appropriate D2P(current_central_density)
      if currentmethod==0 :
          P.append([DensityToPressure2(rho[0][0])])
          P.append([DensityToPressure2(rho[1][0])])
      if currentmethod==1:
          P.append([DensityToPressure3(rho[0][0])])
          P.append([DensityToPressure3(rho[1][0])])

      flags[currentmethod].append([0,0,0])
      for i in range(1,100000): #Arbitrary number of iterations to ensure enough
          if(rho[0][i-1]<=2e18 and rho[0][i-1]>2e17 and boundOne==False):
              flags[currentmethod][-1][0]=(r[0][i-1])/1000
              boundOne=True
          if(rho[0][i-1]<2e17 and rho[0][i-1]>4.3e14 and boundTwo==False):
              flags[currentmethod][-1][1]=(r[0][i-1])/1000
              boundTwo=True
          if(rho[0][i-1]<4.3e14 and rho[0][i-1]>1e9 and boundThree==False):
              flags[currentmethod][-1][2]=(r[0][i-1])/1000
              boundThree=True
          #Use the previous mass,radius and density to calculate the next one

          if(continueSecondStar):
              mResults= rk4([m[0][i-1],m[1][i-1]], Mderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              mResults= rk4([m[0][i-1],m[1][-1]], Mderiv ,[r[0][i-1],r[1][-1]],h,[rho[0][i-1],rho[1][-1]],[m[0][i-1],m[1][-1]])

          """
          #RK5 for comparison
          if(continueSecondStar):
              mResults= rk5([m[0][i-1],m[1][i-1]], Mderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              mResults= rk5([m[0][i-1],m[1][-1]], Mderiv ,[r[0][i-1],r[1][-1]],h,[rho[0][i-1],rho[1][-1]],[m[0][i-1],m[1][-1]])
          """
          """
          #Euler for comparison
          if(continueSecondStar):
              mResults= Euler([m[0][i-1],m[1][i-1]], Mderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              mResults= Euler([m[0][i-1],m[1][-1]], Mderiv ,[r[0][i-1],r[1][-1]],h,[rho[0][i-1],rho[1][-1]],[m[0][i-1],m[1][-1]])
          """
          (r1,m1)=mResults[0]
          (r2,m2)=mResults[1]

          #Use the previous mass,radius,density and pressure to calculate the next one

          if(continueSecondStar):
              pResults = rk4([P[0][i-1],P[1][i-1]], Pderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              pResults = rk4([P[0][i-1],P[1][-1]], Pderiv ,[r[0][-1],r[1][-1]],h,[rho[0][-1],rho[1][-1]],[m[0][i-1],m[1][-1]])

          """
          #RK5 for comparison
          if(continueSecondStar):
              pResults = rk5([P[0][i-1],P[1][i-1]], Pderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              pResults = rk5([P[0][i-1],P[1][-1]], Pderiv ,[r[0][-1],r[1][-1]],h,[rho[0][-1],rho[1][-1]],[m[0][i-1],m[1][-1]])
          """
          """
          #Euler for comparison
          if(continueSecondStar):
              pResults = Euler([P[0][i-1],P[1][i-1]], Pderiv ,[r[0][i-1],r[1][i-1]],h,[rho[0][i-1],rho[1][i-1]],[m[0][i-1],m[1][i-1]])
          else:
              pResults = Euler([P[0][i-1],P[1][-1]], Pderiv ,[r[0][-1],r[1][-1]],h,[rho[0][-1],rho[1][-1]],[m[0][i-1],m[1][-1]])
          """
          (r1,P1)=pResults[0]
          (r2,P2)=pResults[1]

          if(P1<=0):#Boundary condition to know when to stop since edge is reached
              #interpolate the pressure and mass to where P=0
              #where RE is radius error due to interploation (newton)
              (RF1,MF1,RE1) = Extrap(P[0][i-1],P[0][i-2],m[0][i-1],m[0][i-2],r[0][i-1],r[0][i-2],h[0])
              #print("Final Mass for h %f : %2.3e"%(h[0],MF1))
              #print("Final Radius for h %f: %2.4e"%(h[0],RF1))
              #print("Final Interpolation Radius Error for h %f: %f"%(h[0],RE1))
              #Append the final star mass and radius for the given method
              MFE[currentmethod].append((abs(MF1-m[1][-1])/15)/(1.989*10**30))
              RFE[currentmethod].append((pow(pow(abs(RF1/1000-r[1][-1]/1000)/15,2) + pow(RE1,2),0.5)))
              MSolarFinal[currentmethod].append((MF1/(1.989*10**30)))
              RFinal[currentmethod].append(RF1/1000)
              #this will be completed after the rk4 double step.
              #rk4 double step error = (abs(mstep-mdoublestep))/(2^n-1) where n is 4 for the rk4.
              #rk5 and euler methods will be used to compare results.
              #Append the final star mass and radius for the given method
              break

          if(P2<=0 and continueSecondStar):#Boundary condition to know when to stop since edge is reached
              #interpolate the pressure and mass to where P=0
              #where RE is radius error due to interploation (newton)
              (RF2,MF2,RE2) = Extrap(P[1][i-1],P[1][i-2],m[1][i-1],m[1][i-2],r[1][i-1],r[1][i-2],h[1])
              #print("Final Mass for h %f : %2.3e"%(h[1],MF2))
              #print("Final Radius for h %f: %2.4e"%(h[1],RF2))
              #print("Final Interpolation Radius Error for h %f: %f"%(h[1],RE2))
              m[1].append(MF2)
              r[1].append(RF2)
              continueSecondStar=False

          #Append these to be used for the (i+1)th pressure and mass finding
          P[0].append(P1)
          r[0].append(r1)
          m[0].append(m1)

          if(continueSecondStar):
              P[1].append(P2)
              r[1].append(r2)
              m[1].append(m2)

          #Decide which P2D method to use depending on the current method being investigated
          if(rho[0][i-1]>breakpoint and currentmethod==0):
              rho[0].append(PressureToDensity2(P1))
          elif(rho[0][i-1]<=breakpoint and currentmethod==0):
              rho[0].append(PressureToDensity1(P1))
          elif(currentmethod==1):
              rho[0].append(PressureToDensity3(P1))

          if(continueSecondStar):
              if(rho[1][i-1]>breakpoint and currentmethod==0):
                  rho[1].append(PressureToDensity2(P2))
              elif(rho[1][i-1]<=breakpoint and currentmethod==0):
                  rho[1].append(PressureToDensity1(P2))
              elif(currentmethod==1):
                  rho[1].append(PressureToDensity3(P2))
print("\n")
#A graph of Mass vs Radius for the whole family,for all methods
fig = plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    x2=(family*10**15)
    yerror=MFE[i]
    xerror=RFE[i]
    plt.xlabel("Radius/km")
    plt.ylabel("Solar Masses")
    plt.errorbar(x1, y,yerr=yerror,xerr=xerror,fmt="o",label=methodNames[i])
    print("Using method %s the maximum radius is %2.2f km with mass %2.2f M0"%(methodNames[i],max(x1),y[x1.index(max(x1))]))
    print("Using method %s the maximum mass is %2.2f M0 with radius %2.2f km"%(methodNames[i],max(y),x1[y.index(max(y))]))
plt.legend()

#A graph of Mass vs Central Density for the whole family,for all methods
fig2=plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    yerror=MFE[i]
    x2=(family*10**15)
    plt.errorbar(x2, y,yerr=yerror,fmt="o",label=methodNames[i])
    plt.xlabel("Central Density/ kgm^-3")
    plt.ylabel("Solar Masses")

plt.legend()
#A graph of Radius vs Central Density for the whole family,for all methods

fig3=plt.figure()
for i in range(methods):
    x1=RFinal[i]
    y=MSolarFinal[i]
    y2=RFinal[i]
    yerror=MFE[i]
    x2=(family*10**15)
    plt.errorbar(x2, y2,yerr=yerror,fmt="o",label=methodNames[i])
    plt.xlabel("Central Density/ kgm^-3")
    plt.ylabel("Radius/km")
plt.legend()

fig4, ax = plt.subplots()

r1=(flags[0][35][0])
r2=(flags[0][35][1])
r3=(flags[0][35][2])
r4=RFinal[0][35]

ax.annotate('Core',rotation=-45,xy=(r1/2, r1/2))
ax.annotate('Neutron Liquid',xy=(r2/2, r2/2),rotation=-45)
ax.annotate('Inner Crust',xy=(2/3 *r3, r3/2),rotation=-45)
ax.annotate("Outer Crust",
            xy=(r4*math.cos(math.radians(45)),r4*math.sin(math.radians(45))), xycoords='data',
            xytext=(r4-3, r4-1), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax.add_patch(Wedge((0,0), r4,0,90,color="green"))
ax.add_patch(Wedge((0,0), r3,0,90,color="yellow"))
ax.add_patch(Wedge((0,0), r2,0,90,color="orange"))
ax.add_patch(Wedge((0,0), r1,0,90,color="red"))
ax.set_title("Method: %s with Central Density: %2.3e"%(methodNames[0],family[35]*1e17))
ax.plot()
plt.xlabel("Radius/km")
plt.ylabel("Radius/km")

fig5, ax2 = plt.subplots()

r1=(flags[1][35][0])
r2=(flags[1][35][1])
r3=(flags[1][35][2])
r4=RFinal[1][35]

ax2.annotate('Core',xy=(r1/2, r1/2),rotation=-45)
ax2.annotate('Neutron Liquid',xy=(r2/2, r2/2),rotation=-45)
ax2.annotate("Outer Crust",
            xy=(r4*math.cos(math.radians(45)),r4*math.sin(math.radians(45))), xycoords='data',
            xytext=(r4-3, r4-1), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax2.add_patch(Wedge((0,0), r4,0,90,color="green"))
ax2.add_patch(Wedge((0,0), r3,0,90,color="yellow"))
ax2.add_patch(Wedge((0,0), r2,0,90,color="orange"))
ax2.add_patch(Wedge((0,0), r1,0,90,color="red"))
ax2.plot()
ax2.set_title("Method: %s with Central Density: %2.3e"%(methodNames[1],family[35]*1e17))
plt.xlabel("Radius/km")
plt.ylabel("Radius/km")
plt.show()
