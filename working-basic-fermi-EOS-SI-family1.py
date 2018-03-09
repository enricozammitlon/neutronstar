import numpy as np
import matplotlib.pyplot as plt
c=299792458
G = 6.67*10**-11
hbar = 1.055*10**(-34)
MNEUTRON = 1.67*10**(-27)

MSolarFinal =[]
RFinal =[]

family = np.arange((3.7), (3.7*10**2),20)
for rhoIn in family:

    h = 100.0
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
    '''
    print(DensityToPressure1(rho[0]))
    print(PressureToDensity1(DensityToPressure1(rho[0])))
    
    print(DensityToPressure2(rho[0]))
    print(PressureToDensity2(DensityToPressure2(rho[0])))
    '''
    def Pderiv(P,r,rho,m):
        
        return (-G*m*rho/pow(r,2))
    
    def Mderiv(M,r,rho,b):
        return (4*rho*np.pi*pow(r,2))
    
    P.append(DensityToPressure2(rho[0]))
    #P.append(500)
    for i in range(1,10000000):
        if i < 5:
            print("Mass:    \t%2.6e"%(m[i-1]))
            print("Pressure:\t%2.6e"%(P[i-1]))
            print("New Density:\t%e"%(rho[i-1]))
            print("New R :\t%f"%(r[i-1]))
            (r1,m1)= rk4(m[i-1], Mderiv ,r[i-1],h,rho[i-1],1)
            (r1,P1) = rk4(P[i-1], Pderiv ,r[i-1],h,rho[i-1],m[i-1])
            if(P1<=0):
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
        
    
    '''   
    fig = plt.figure()
    x=r
    y1=MSOL
    y2=P
    plt.xlabel("radius meters")
    plt.ylabel("SOLAR MASSES")
    
    plt.plot(x,y1)
    
    fig2=plt.figure()
    plt.plot(x,y2)
    plt.xlabel("radius meters")
    plt.ylabel("pressure Nm^-2")
    plt.show()
    '''
print (MSolarFinal, RFinal)

fig = plt.figure()
x1=RFinal
y=MSolarFinal
x2=(family*10**15)
plt.xlabel("radius KM")
plt.ylabel("SOLAR MASSES")

plt.plot(x1,y)

fig2=plt.figure()
plt.plot(x2,y)
plt.xlabel("central density")
plt.ylabel("SOLAR MASSES")
plt.show()

    
    
    
    
    
    
    