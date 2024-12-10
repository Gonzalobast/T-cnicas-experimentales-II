from numpy import *
import numpy as np
from numpy.linalg import *

import scipy.optimize as so

import matplotlib.pyplot as plt

xmax=np.array([0.2855694407477436, 0.8862474407477436, 1.4869254407477437, 2.0876034407477437, 2.688281440747744, 3.2889594407477443])

xmin= np.array([0.0, 0.5711388814954872, 1.1422777629909744, 1.7134166444864616, 2.2845555259819488, 2.855694407477436])

xt = np.array([0.0, 0.2855694407477436, 0.5711388814954872, 0.8862474407477436, 1.1422777629909744, 1.4869254407477437, 1.7134166444864616, 2.0876034407477437, 2.2845555259819488, 2.688281440747744, 2.855694407477436, 3.2889594407477443])


y6=np.array([-18,12.4,-9.2,6.4,-4.8,3.4,-2.2,1.6,-0.9,0.6,-0.4,0.2])

y6max = np.array([12.4, 6.4, 3.4, 1.6, 0.6, 0.2])
y6min = np.array([-18.0, -9.2, -4.8, -2.2, -0.9, -0.4])

yerror = 0.2

#Definimos la función con los parámetros a ajustar

def fun(t,a,b): 
    y = a*exp(-b*t)
    return y

def theta(t,b,c,d):
    y =-18*np.exp(-b*t)*np.cos(c*t+d)
    return y

# Introducimos valores iniciales arbitrarios para los parámetros que hay que ajustar
# par = [w,gamma,A,phi,y0]

par = [0.1,1]
par2   = [0.1,1,1]

solmax = so.curve_fit(fun,xmax,y6max,p0=(par))
solmin = so.curve_fit(fun,xmin,y6min,p0=(par))
sol2 = so.curve_fit(theta,xt,y6,p0=(par2))

# La salida tiene dos arrays: sol[0] con los parámetros  
amax,bmax = solmax[0] #valores estimados de los parametros
amin, bmin = solmin[0]
b2,c,d = sol2[0]

# y sol[1] con las incertidumbres. Ésta es una matriz de covarianza. Sí solo queremos las varianzas tomamos la diagonal
samax,sbmax = np.sqrt(np.diag(solmax[1]))
samin,sbmin = np.sqrt(np.diag(solmin[1]))
sb2,sc,sd = np.sqrt(np.diag(sol2[1]))
#Incertidumbres de los parámetros

# la función estimada
xprueba = np.arange(0,4,0.01)
yEstmax = fun(xprueba,amax,bmax)
yEstmin = fun(xprueba,amin,bmin)
yest2 = theta(xprueba,b2,c,d)

print(f'Para los maximos el valor de a es:{amax} y su incertidumbre:{samax} y el de b es: {bmax} y su incertidumbre:{sbmax}')
print(f'Para los minimos el valor de a es:{amin} y su incertidumbre:{samin} y el de b es: {bmin} y su incertidumbre:{sbmin}')
print(f'b={b2} inc:{sb2}, c={c} inc:{sc}, d={d} inc:{sd}')


#plt.plot(xprueba,yEst,'-',color='royalblue',label='Ajuste a una función $\u03B8_{max}=ae^{-bt}$')
plt.plot(xprueba,yest2,'-',color='red',label='Ajuste a una función $\u03B8=\u03B8_{0}e^{-bt} cos(ct+d)$')
plt.plot(xprueba,yEstmax,'-',color='royalblue',label='Ajuste a una función $\u03B8_{max}=ae^{-bt}$')
plt.plot(xprueba,yEstmin,'-',color='royalblue')
plt.errorbar(xmax,y6max,yerror,marker='o',linestyle='None',color='orange',label='Puntos experimentales')
plt.errorbar(xmin,y6min,yerror,marker='o',linestyle='None',color='orange')
plt.xlabel('t(s)')
plt.ylabel('$\u03B8$')
#plt.errorbar(x,y,yerror,marker='.',linestyle='None',color='orange',label='Puntos experimentales')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
