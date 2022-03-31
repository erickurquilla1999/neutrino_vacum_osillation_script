import numpy as np
import matplotlib.pyplot as plt

rho_electron=[]
time=[]

for i in range(0,20000):
    data = np.loadtxt("output/step_"+str(i)+".dat", unpack = True)
    rho_electron.append(data[0][0])
    time.append(data[2][0])

theta_12=33.82*np.pi/180
theta_13=8.61*np.pi/180
theta_23=48.3*np.pi/180
E=10.0e6

s12=np.sin(theta_12)
c12=np.cos(theta_12)
s13=np.sin(theta_13)
c13=np.cos(theta_13)
s23=np.sin(theta_23)
c23=np.cos(theta_23)
hbar=1.054571817e-34

phi32=2.449e-3/(4*E)
phi31=2.449e-3/(4*E)
phi21=7.39e-5/(4*E)

#print(phi32)
#print(phi31)
#print(phi21)
#print(np.sin(2*theta_13)**2*c12**2)
#print(np.sin(2*theta_13)**2*s12**2)
#print(c13**4*np.sin(2*theta_12)**2)

Pee=[]
ti=[]

for t in np.linspace(0,(1/1.8475e-12)*2*np.pi*0.6,10000):
    ti.append(t)
    Pee.append(1-np.sin(2*theta_13)**2*c12**2*np.sin(phi31*t)**2-np.sin(2*theta_13)**2*s12**2*(np.sin(phi32*t))**2-(c13**4)*np.sin(2*theta_12)**2*(np.sin(phi21*t))**2)

plt.plot(ti,Pee,label=r'Analitic solution')
plt.plot(time,rho_electron,label=r'Numerical solution')
plt.ylabel("Probability")
plt.xlabel("time")
plt.legend()
plt.show()
