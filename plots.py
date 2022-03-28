import numpy as np

rho_electron=[]
time=[]

for i in range(0,100):
    data = np.loadtxt("step_"+str(i)+".dat", unpack = True)
    rho_electron.append(data[0][0])
    time.append(data[2][0])

import matplotlib.pyplot as plt
plt.plot(time,rho_electron)
plt.show()

M_PI=np.pi
#m1=0.049*1.782662e-36
m1=0.049
m2=0*1.782662e-36
m3=0*1.782662e-36
theta_12=1e-6*M_PI/180
theta_13=48.3*M_PI/180
theta_23=8.61*M_PI/180
#E=10.0e6*1.60218e-19
E=10.0e6
#c=299792458
c=1
s12=np.sin(theta_12)
c12=np.cos(theta_12)
s13=np.sin(theta_13)
c13=np.cos(theta_13)
s23=np.sin(theta_23)
c23=np.cos(theta_23)
hbar=1.054571817e-34
phi32=(m3-m2)**2*c**4/(4*E)
phi31=(m3-m1)**2*c**4/(4*E)
phi21=(m2-m1)**2*c**4/(4*E)

print(phi32)
print(phi31)
print(phi21)


Pee=[]
ti=[]
for t in np.linspace(0,(1/6.002500000000001e-11)*6,1000):
    ti.append(t)
    Pee.append(1-np.sin(2*theta_13)**2*(c12**2*np.sin(phi31*t)**2+s12**2*np.sin(phi32*t)**2)-c13**4*np.sin(2*theta_12)**2*np.sin(phi21*t)**2)

plt.plot(ti,Pee)
plt.show()
