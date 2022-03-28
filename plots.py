import numpy as np

rho_electron=[]
time=[]

for i in range(0,100):
    data = np.loadtxt("/home/centroescolarjuanabarrera/tesis/v_o_s/output/step_"+str(i)+".dat", unpack = True)
    rho_electron.append(data[0][0])
    time.append(data[2][0])

import matplotlib.pyplot as plt
plt.plot(time,rho_electron)
plt.show()

