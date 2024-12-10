from numpy import *
import numpy as np
from numpy.linalg import *

import scipy.optimize as so

import matplotlib.pyplot as plt


from IPython.display import display, Latex

import pandas as pd

x = np.array([0.0, 0.300339, 0.600678, 0.9010170000000001, 1.201356, 1.5016950000000002, 1.8020340000000001, 2.102373, 2.402712, 2.7030510000000003, 3.0033900000000004, 3.303729, 3.6040680000000003, 3.9044070000000004])
nx = len(x)

y = np.array([18,14.8,12.4,10.2,8.6,7,5.8,4.6,3.8,3.2,2.8,2.2,1.8,1.2])
ny = len(y)

#the estimated type B error un y

yerror = 0.2

#Definimos la función con los parámetros a ajustar

def fun(t,a,b): 
    y = a*exp(-b*t)
    return y

# Introducimos valores iniciales arbitrarios para los parámetros que hay que ajustar
# par = [w,gamma,A,phi,y0]

par   = [0.1,1]

sol = so.curve_fit(fun,x,y,p0=(par))

# La salida tiene dos arrays: sol[0] con los parámetros  
a,b = sol[0] #valores estimados de los parametros

# y sol[1] con las incertidumbres. Ésta es una matriz de covarianza. Sí solo queremos las varianzas tomamos la diagonal
sa,sb = np.sqrt(np.diag(sol[1])) #Incertidumbres de los parámetros

# la función estimada
xprueba = np.arange(0,8,0.01)
yEst=fun(xprueba,a,b)

print(f'El valor de a es:{a} y su incertidumbre:{sa} y el de b es: {b} y su incertidumbre:{sb}')

curva = fun(xprueba,a,b)


plt.plot(xprueba,yEst,'-',color='royalblue',label='Ajuste a una función $\u03B8_{max}=ae^{-bt}$')
plt.xlabel('t(s)')
plt.ylabel('$\u03B8_{max}$')
#plt.plot(xprueba,curva,'-',color='royalblue')
plt.errorbar(x,y,yerror,marker='.',linestyle='None',color='orange',label='Puntos experimentales')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
