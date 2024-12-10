import numpy as np
import matplotlib.pyplot as plt
from math import *
import scipy.optimize as so

from fig_config import (
    add_grid,
    figure_features
)

figure_features()

def regresionSimple(x,y):
    """Axusta os datos dos vectore x e y a unha resta dada pola ec. y=a + bx
    Parametros:
    x vector con medidas da magnitud x
    y vector con medidas da magnitud y
    Devolve:
    a coeficiente a
    b coeficiente b
    sa incerteza de a
    sb incerteza de b
    r coeficiente de regresion lineal """
    n=len(x)
    sx=sum(x); sy=sum(y); xx=np.dot(x,x); yy=np.dot(y,y); xy=np.dot(x,y)
    denom=(n*xx - sx**2)
    b=(n*xy - sx*sy)/denom
    a=(xx*sy - sx*xy)/denom
    s=sqrt(sum((y-a-b*x)**2)/(n-2))
    sa=s*sqrt(xx/(n*xx-sx**2))
    sb=s*sqrt(n/(n*xx-sx**2))
    r=(n*xy-sx*sy)/sqrt((n*xx-sx**2)*(n*yy-sy**2))
    return [a,b, sa, sb, r, s]

### movimiento de nutacion

## (t,omega), omega en medias vueltas por segundo, t es el tiempo de 4 cabeceos
serie1n = np.array([[3.85,4.38,5.41,5.67,6.16],[15.9,12.9,10.9,9.6,8.3]])
serie2n = np.array([[5.03,5.44,6.60,6.93,7.77],[9.9,9.1,7.7,7.1,6.4]])
serie3n = np.array([[3.64,4.06,4.50,4.80,5.01],[15.3,14.2,12.4,11.4,10.8]])

def omega_nut(t):
    st = 0.5
    n = 4
    omega_n = 2*pi*n/t
    s_omega_n = 2*pi*n*st/(t**2)
    return [omega_n,s_omega_n]

def omega(f):
    sf = 0.5
    omega = pi*f #pasamos a vueltas completas ya
    s_omega=pi*sf
    return [omega,s_omega]

def media_ponderada(x,sx):
    w = 1/sx**2
    media = np.dot(w,x)/sum(w)
    inc_media = 1/sqrt(sum(w))
    return [media,inc_media]

r = (25.5/2)/100 #r=diam/2
sr = 0.05/100
M= 1.741

I3_teor = M*r**2/2
s_I3_teor = M*r*sr

def i_nut(b,sb):
    i1 = I3_teor/b
    s_i1 = sqrt((s_I3_teor/b)**2+(I3_teor*sb/b**2)**2)
    return [i1,s_i1]

#Regresion 1

omega_nut1 = omega_nut(serie1n[0])
omega1 = omega(serie1n[1])

reg_nut1 = regresionSimple(omega1[0],omega_nut1[0])
i1_nut = i_nut(reg_nut1[1],reg_nut1[3])

#Regresion 2

omega_nut2 = omega_nut(serie2n[0])
omega2 = omega(serie2n[1])

reg_nut2 = regresionSimple(omega2[0],omega_nut2[0])
i2_nut = i_nut(reg_nut2[1],reg_nut2[3])

#Regresion 3

omega_nut3 = omega_nut(serie3n[0])
omega3 = omega(serie3n[1])

reg_nut3 = regresionSimple(omega3[0],omega_nut3[0])
i3_nut = i_nut(reg_nut3[1],reg_nut3[3])

momentos = np.array([i1_nut[0],i2_nut[0],i3_nut[0]])
s_momentos = np.array([i1_nut[1],i2_nut[1],i3_nut[1]])

i_nut = media_ponderada(momentos,s_momentos)

print(i_nut)



fig = plt.figure(figsize=(8,6))
ax = plt.axes(xlim=(15,55),ylim=(0,9))

ax.set_xlabel(r'$\omega(rad\cdot s^{-1})$')
ax.set_ylabel(r'$\Omega_{n}(rad\cdot s^{-1})$')

add_grid(ax,locations=(0.5,1,1,2))

ax.set_xticks(np.arange(15,55.1,5))
ax.set_yticks(np.arange(0,9.1,1))

ax.errorbar(omega1[0],omega_nut1[0],xerr=pi/2,yerr=omega_nut1[1],linestyle='None',marker='.',markersize=10,color='b',label='Serie 1')
ax.plot(omega1[0],reg_nut1[0]+reg_nut1[1]*omega1[0],'--',color='b')
ax.errorbar(omega2[0],omega_nut2[0],xerr=pi/2,yerr=omega_nut2[1],linestyle='None',marker='.',markersize=10,color='g',label='Serie 2')
ax.plot(omega2[0],reg_nut2[0]+reg_nut2[1]*omega2[0],'--',color='g')
ax.errorbar(omega3[0],omega_nut3[0],xerr=pi/2,yerr=omega_nut3[1],linestyle='None',marker='.',markersize=10,color='r',label='Serie 3')
ax.plot(omega3[0],reg_nut3[0]+reg_nut3[1]*omega3[0],'--',color='r')

ax.legend(loc='lower right')

plt.show()