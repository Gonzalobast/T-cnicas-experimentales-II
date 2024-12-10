import numpy as np
from math import *
import matplotlib.pyplot as plt
from numpy.linalg import *
import scipy.optimize as so


frecuencias = np.array([[14.2, 28.5, 43.5, 57.5, 71.4, 83.6],
    [20.1, 40.2, 60.3, 80.2, 100, 120],
    [31.7, 63.6, 95, 126.7, 158.9, 190.7],
    [37.4, 74.7, 112, 149.4, 187.3, 224.2],
    [44.5, 89, 133.4, 177.9, 222.4, 266.9]])

d = np.array([1,675, 0.838, 0.558, 0.419, 0.335,0.279])
landa = np.array([3.350, 1.675, 1.117, 0.838, 0.670,0.558])

T = np.array([[0.070422535, 0.035087719, 0.022988506, 0.017391304, 0.014005602, 0.011961722],
    [0.049751244, 0.024875622, 0.016583748, 0.012468828, 0.01, 0.008333333],
    [0.031545741, 0.01572327, 0.010526316, 0.00789266, 0.006293266, 0.005243838],
    [0.026737968, 0.013386881, 0.008928571, 0.00669344, 0.005339028, 0.004460303],
    [0.02247191, 0.011235955, 0.007496252, 0.005621135, 0.004496403, 0.003746722]])

def regresionSimpleSenTermoIndependente(x,y):
    """Axusta os datos dos vectore x e y a unha resta dada pola ec, y= bx
    Parametros:
    x vector con medidas da magnitud x
    y vector con medidas da magnitud y
    Devolve:
    b coeficiente b
    sb incerteza de b
    r coeficiente de regresion lineal """
    n=len(x)
    xy=np.dot(x,y); xx=np.dot(x,x); xy=np.dot(x,y); yy=np.dot(y,y)
    b=xy/xx
    s=sqrt(sum((y-b*x)**2)/(n-1))
    sb=s/sqrt(xx)
    r=xy/sqrt(xx*yy)
    """print('b = ',b)
    print('s(b) =',sb)
    print('r = ',r)
    print('s = ',s)"""
    return [b, sb, r, s]

masa1 = regresionSimpleSenTermoIndependente(T[0],landa)
masa2 = regresionSimpleSenTermoIndependente(T[1],landa)
masa3 = regresionSimpleSenTermoIndependente(T[2],landa)
masa4 = regresionSimpleSenTermoIndependente(T[3],landa)
masa5 = regresionSimpleSenTermoIndependente(T[4],landa)

x1 = np.linspace(0,T[0,0],500)
x2 = np.linspace(0,T[1,0],500)
x3 = np.linspace(0,T[2,0],500)
x4 = np.linspace(0,T[3,0],500)
x5 = np.linspace(0,T[4,0],500)

"""plt.plot(T[0],landa,'o',color='b')
plt.plot(x1,masa1[0]*x1,color='b',label=r'$m_{1}$')
plt.plot(T[1],landa,'o',color='g')
plt.plot(x2,masa2[0]*x2,color='g',label=r'$m_{2}$')
plt.plot(T[2],landa,'o',color='r')
plt.plot(x3,masa3[0]*x3,color='r',label=r'$m_{3}$')
plt.plot(T[3],landa,'o',color='c')
plt.plot(x4,masa4[0]*x4,color='c',label=r'$m_{4}$')
plt.plot(T[4],landa,'o',color='m')
plt.plot(x5,masa5[0]*x5,color='m',label=r'$m_{5}$')
plt.xlabel('T(s)')
plt.ylabel(r'$\lambda$(m)')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()"""


#Incertidumbre del periodo

s_f = 0.1

def incT(x):
    incT = s_f/(x**2)
    return incT

incT = incT(frecuencias)
nx,ny = incT.shape

for i in range(nx):
    for j in range(ny):
        incT[i,j] ="{:.7f}".format(incT[i,j])

### Densidad lineal

#Calculamos c^2

c = np.array([masa1[0],masa2[0],masa3[0],masa4[0],masa5[0]])
s_c = np.array([masa1[1],masa2[1],masa3[1],masa4[1],masa5[1]])

c_cuadrado = c**2
s_c_cuadrado = s_c**2

masas = np.array([100,198,504,704,1002])/1000

tau = masas*9.8
s_tau = 9.8 #Incertidumbre 

reg_dens_lineal = regresionSimpleSenTermoIndependente(tau,c_cuadrado)

"""plt.plot(tau,c_cuadrado,'o',color='orange',label='Datos experimentales')
plt.plot(tau,tau*reg_dens_lineal[0],color='black',linestyle='--',label='Recta de regresión')
plt.xlabel(r'$\tau(N)$')
plt.ylabel(r'$c^{2}(m^{2}\cdot s^{-2})$')
plt.grid(True)
plt.legend(loc='upper left')
plt.show()"""

rho_l = 1/reg_dens_lineal[0]
s_rho_l = reg_dens_lineal[1]/(reg_dens_lineal[0]**2)

rho_v = 7850

d = 2*sqrt(rho_l/(pi*rho_v))
s_d = s_rho_l/(pi*rho_v)*1/(sqrt(rho_l/(pi*rho_v)))

print(f'd: {d} y s(d): {s_d}')