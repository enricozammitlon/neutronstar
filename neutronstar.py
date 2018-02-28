import math
import numpy as np

global G
global MSOLAR
global c
c=299792458
G=197.327*6.67259*(10**-45)
MSOLAR=1.1157467*(10**60)
rhoS=[((3.7*10**17)*10**(-45)*(3*10**8)**2)/(1.6*10**(-13))]

def initConstants(rhoS):
    rZero=[]
    mZero=[]
    check=[]
    for rhos in rhoS:
        rZero.append((10**-19)*(1/(np.sqrt(rhos*G*4*np.pi))))
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))
        check.append((G*4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2)))/(1/(np.sqrt(rhos*G*4*np.pi))))
    return rZero,mZero,check
print(initConstants(rhoS))

#check caclulation should without MSOLAR or Km correction = 1
