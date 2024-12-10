import numpy as np
import matplotlib.pyplot as plt

## rho(x)=0.02026x^2-0.0435x+0.8341

rho = np.array([1])

x = np.array([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0])
rho29 = np.array([0.9896,0.9801,0.9663,0.9373,0.9192,0.8982,0.8742,0.8546,0.8429,0.8326])

"""pol = np.polyfit(x,rho29,2)

print(pol)"""

"""a = 0.00261568
b = 0.605354004
c = 0.8324713564"""

a=-28.048
b=55.761
c=-26.75



T_eb = np.array([99.8,95.2,92.5,90.4,88.7,87.4,77.9,78.3,78.9,79.7,80.2,80.7,81.3,81.8,82.3,83.6,84.6,85.5])

#x(rho)

def x_rho(rho): 
    x1 = (-b+np.sqrt((b**2)-4*a*(c-rho)))/(2*a)
    x2 = (-b-np.sqrt((b**2)-4*a*(c-rho)))/(2*a)
    print(x1)
    return x1

rho30_liq = np.array([0.9879,0.9357,0.914,0.9051,0.8925,0.8802,0.8093,0.8155,0.8202,0.827,0.8323,0.8395,0.8447,0.8515,0.8554,0.861,0.8687,0.8717])
rho30_vap = np.array([0.9945,0.9876,0.9806,0.9774,0.9688,0.9662,0.8039,0.822,0.8395,0.8602,0.8752,0.8816,0.9002,0.9142,0.925,0.9448,0.9542,0.9606])
rho28_liq = np.array([0.9885,0.9371,0.9174,0.9128,0.8978])
rho28_vap = np.array([0.9952,0.9883,0.9815,0.9781,0.9673])

xprueba30_liq = a*rho30_liq**2+b*rho30_liq+c
print(xprueba30_liq)

x30_liq = np.array(x_rho(rho30_liq))
x28_liq = x_rho(rho28_liq)
x30_vap = np.array(x_rho(rho30_vap))
x28_vap = x_rho(rho28_vap)

y1= a*rho30_liq**2+b*rho30_liq+c
y2 = a*rho30_vap**2+b*rho30_vap+c


#Representación para 30ºC

plt.plot(y1,T_eb,'o')
plt.plot(y2,T_eb,'o')

plt.show()





