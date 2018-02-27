import math
import numpy as np

global G
global MSOLAR
global c
c=299792458
G=197.327*6.67259*(10**-45)
MSOLAR=1.1157467*(10**60)
rhoS=[(32.7*10**17)*1.6*10/(9*10**16)]

def initConstants(rhoS):
    rZero=[]
    mZero=[]
    for rhos in rhoS:
        rZero.append(((10)**(-19))/(np.sqrt(rhos*G*4*np.pi)))
        mZero.append((4*np.pi*rhos)/(((rhos*G*4*np.pi)**(3/2))*MSOLAR))
    return rZero,mZero

print(initConstants(rhoS))
