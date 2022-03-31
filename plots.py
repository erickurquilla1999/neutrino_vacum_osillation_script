import numpy as np
import matplotlib.pyplot as plt

#simulation
rho_electron=[]
time=[]

for i in range(0,10000):
    data = np.loadtxt("output/step_"+str(i)+".dat", unpack = True)
    rho_electron.append(data[0][0])
    time.append(data[2][0])

#constant

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

m1=1.0
m2=10.0
m3=10.1

phi32=(m3**2-m2**2)/(4*E)
phi31=(m3**2-m1**2)/(4*E)
phi21=(m2**2-m1**2)/(4*E)

#analityc data

Pee=[]
ti=[]

for t in np.linspace(0,(1/5.024999999999977e-08)*2*np.pi*0.6,10000):
    ti.append(t)
    Pee.append(1-np.sin(2*theta_13)**2*c12**2*np.sin(phi31*t)**2-np.sin(2*theta_13)**2*s12**2*(np.sin(phi32*t))**2-(c13**4)*np.sin(2*theta_12)**2*(np.sin(phi21*t))**2)

#making plots

fig,ax=plt.subplots()
plt.plot(ti,Pee,label=r'Analitic solution')
plt.plot(time,rho_electron,label=r'Numerical solution')
plt.ylabel(r"$\rho_{ee}$")
plt.xlabel(r"Time (natural units)")
plt.xlim([-0.1e7, 3e7])
plt.legend()
plt.savefig("vaccum_neutrino_flavor_oscilattion_plot.pdf")
plt.show()
