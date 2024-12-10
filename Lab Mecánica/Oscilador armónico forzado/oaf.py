import numpy as np
from math import *
from numpy.linalg import *

import scipy.optimize as so

import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

def VariasMedidas(x):   
    k = 2 #Factor de cobertura
    
    sb = 0.5
    
    media1 = np.mean(x)
    desviacion1 = np.std(x)
    n = len(x)
    
    #Vamos a crear el intervalo de confianza
    
    limInferior = media1 - k*desviacion1
    limSuperior = media1 + k*desviacion1
    
    datos=[];excluidos=[]
    
    for i in range(n):
        if ((x[i] >= limInferior) and (x[i] <= limSuperior)):
            datos.append(x[i])  #Lista con los datos filtrados
        else:
            excluidos.append(x[i])  #Lista con los datos excluidos
            
    media2 = np.mean(datos)    #Nuevo valor de la media
    
    sa = np.std(datos)/sqrt(len(datos)) #Incertidumbre de tipo a de la media
    
    sc = sqrt(sa**2 + sb**2)    #Incertidumbre combinada

    """print(f'Datos de la muestra:')
    print(f'media inicial: {media1}')
    print(f'Desv tipica inicial: {desviacion1}')
    print(f'sb: {sb}')
    print(f'Datos excluidos: {excluidos}')
    print(f'Datos finales: {datos}')
    print(f'Media final: {media2}')
    print(f'Desv tipica final de la muestra: {np.std(datos)}')
    print(f'sa de la media {sa}')
    print(f'sc de la muestra: {sc}')"""
    
    return (sa,sb,sc,media2,np.std(datos),excluidos,datos,media1,desviacion1)

def frecuencia(t,st,n):
    omega1 = 2*pi*n/t
    T = 2*pi/omega1
    s_omega1 = 2*pi*n*st/(t**2)
    print('Para t=',t,' la frecuencia es ',omega1,' y la incertidumbre es ',s_omega1,'El periodo es: ',T)
    return(omega1,s_omega1,T)

#Valores del tiempo medidos para calcular el periodo y omega

I0 = np.array([56.73, 56.83, 56.22, 56.67])
I3 = np.array([18.7,18.6,18.6,18.76,18.80])
I6 = np.array([12.72,12.67,12.21,12.86,12.34])
I9 = np.array([7.41,7.49,7.48,7.61,7.47])
I11 = np.array([5.52,5.44,5.65,5.64,5.71])

n0,n3,n6,n9,n11 = 30,10,7,4,3

def media(I):
    medias = VariasMedidas(I)
    return (medias[3],medias[2])

t0 = media(I0)
omega0 = frecuencia(t0[0],t0[1],30)
t3 = media(I3)
omega3 = frecuencia(t3[0],t3[1],10)
t6 = media(I6)
omega6 = frecuencia(t6[0],t6[1],7)
t9 = media(I9)
omega9 = frecuencia(t9[0],t9[1],4)
t11 = media(I11)
omega11 = frecuencia(t11[0],t11[1],3)

omega1 = np.array([omega0[0],omega3[0],omega6[0],omega9[0],omega11[0]])
s_omega1 = np.array([omega0[1],omega3[1],omega6[1],omega9[1],omega11[1]])

#Amplitudes medidas

y0 = np.array([18,17.8, 17, 16.6, 16.2, 15.8, 15.4, 15.2, 14.8, 14.4])
y3 = np.array([18,14.8,12.4,10.2,8.6,7,5.8,4.6,3.8,3.2,2.8,2.2,1.8,1.2])
y6 = np.array([-18,12.4,-9.2,6.4,-4.8,3.4,-2.2,1.6,-0.9,0.6,-0.4,0.2])
y9 = np.array([-18,10.2,-5.2,2.8,-1.2,0.6,-0.2,0.1])
y11 = np.array([-18,8.8,-4,1.8,-0.8,0.4,-0.2])

#Función para obtener los valores de x

def intervalo_t(n,intervalo,T):
    x = []
    for i in range(n):
        x.append(T*i*intervalo)
    return x

x0 = np.array(intervalo_t(len(y0),2,omega0[2]))
x3 = np.array(intervalo_t(len(y3),1,omega3[2]))
x6 = np.array(intervalo_t(len(y6),1/2,omega6[2]))
x9 = np.array(intervalo_t(len(y9),1/2,omega9[2]))
x11 = np.array(intervalo_t(len(y11),1/2,omega11[2]))

x6min = np.array(intervalo_t((len(y6)//2),1,omega6[2]))
x6max = x6min + omega6[2]/2
y6max = y6[y6>0]
y6min = y6[y6<0]

x9min = np.array(intervalo_t((len(y9)//2),1,omega9[2]))
x9max = x9min + omega9[2]/2
y9max = y9[y9>0]
y9min = y9[y9<0]

x11min = np.array(intervalo_t((len(y11)//2),1,omega11[2]))
x11max = x11min + omega11[2]/2
x11min = np.array([0.   , 1.864, 3.728,5.592])

y11max = y11[y11>0]
y11min = y11[y11<0]

def exponencial(x,a,b):
    y = a*e**(-b*x)
    return y

def exp_cos(x,b,c,d):
    y = -18 * e**(-b*x)*np.cos(c*x+d)
    return y

theta_error = 0.2 #Error del ángulo

#Definimos la función que ajuste los datos

def ajuste(fun,par,x,y):
    sol = so.curve_fit(fun,x,y,p0=(par)) #Array de soluciones y incertidumbres
    valores = sol[0]
    inc = np.sqrt(np.diag(sol[1]))
    print(f'Los valores de los parametros a ajustar son {valores}')
    print(f'Las incertidumbres son: {inc}')
    return valores

par0 = [0.1,0.1]
par1=[0.1,1,1]

ajuste11max = ajuste(exponencial,par0,x11max,y11max)
ajuste11min = ajuste(exponencial,par0,x11min,y11min)
fperiodica = ajuste(exp_cos,par1,x11,y11)
"""
xprueba = np.arange(0,10,0.01)
ymax = exponencial(xprueba,ajuste6max[0],ajuste6max[1])
ymin = exponencial(xprueba,ajuste6min[0],ajuste6min[1])
theta = exp_cos(xprueba,fperiodica[0],fperiodica[1],fperiodica[2])


plt.plot(xprueba,ymax,'-',color='royalblue',label='Ajuste a una función $ae^{-bt} $')
plt.plot(xprueba,ymin,'-',color='royalblue')
plt.plot(xprueba,theta,'-',color='red',label='Ajuste a una función $\u03B8_{0}e^{-bt}cos(ct+d)$')
plt.errorbar(x6,y6,theta_error,marker='o',linestyle='None',color='orange',label='Puntos experimentales')
plt.xlabel('t(s)')
plt.ylabel('$\u03B8$')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()"""

#####Graficas de omega0^2

#Array de gamma

f0cuadrado = omega1[0]**2
sf0 = s_omega1[0]*2*omega1[0] #Incertidumbre de omega0

I = np.array([0,0.3,0.6,0.9,1.1])

gamma = np.array([0.006642,0.10221,0.370,0.663,0.8077])
s_gamma = np.array([0.00022,0.00092,0.012,0.16,0.056])

patron = np.array([omega1[0]**2,omega1[0]**2,omega1[0]**2,omega1[0]**2,omega1[0]**2])
s_patron = np.array([2*omega1[0]*s_omega1[0],2*omega1[0]*s_omega1[0],2*omega1[0]*s_omega1[0],2*omega1[0]*s_omega1[0],2*omega1[0]*s_omega1[0]])

comparacion = np.array(omega1[0]**2-(gamma**2)) #Termino omega0^2-gamma^2
s_comparacion = []

for i in range(len(comparacion)):
    inc_i = sqrt((2*omega1[0]*s_omega1[0])**2+(2*gamma[i]*s_gamma[i])**2)
    s_comparacion.append(inc_i)

"""plt.errorbar(I,omega1**2,s_omega1,label=r'$\omega^{2}_{1}$',linestyle='None',marker='s',capsize=5,color='black')
plt.errorbar(I,comparacion,s_comparacion,label=r'$\omega^{2}_{0}-\gamma^{2}$',linestyle='None',marker='s',capsize=5,color='grey')
plt.legend(loc='lower left')
plt.xlabel('I(A)')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
"""

### Gamma(I)

#Hay que ajustar a un polinomio de grado 2



coeff = np.polyfit(I,gamma,2)
gamma_ajuste = np.poly1d(coeff)

recta = []


#Calculo de I_c

x = Symbol('x')

I_c = np.array(solve(coeff[0]*x**2+coeff[1]*x+coeff[2]-omega1[0]))
I_c = I_c[I_c>0]

print(f'Los coeficientes del polinomio son: {coeff}')

sI_c = s_omega1[0]/sqrt(coeff[1]**2+4*coeff[0]*omega1[0]-4*coeff[0]*coeff[2])

print(f'La intensidad crítica es {I_c} y la inc: {sI_c}')

Iprueba = np.linspace(0,I_c,1000)
Iprueba2 = np.linspace(0,3,1000)

for i in range(len(Iprueba)):
    recta.append(omega1[0])

xprueba = []
recta2 = np.linspace(0,omega1[0],1000)

for i in range(len(recta2)):
    xprueba.append(I_c)

print(omega1[0])

####

plt.plot(Iprueba,gamma_ajuste(Iprueba),label=r'$\gamma=aI^{2}+bI+c$')
plt.plot(I,gamma,'o',label='Medidas experimentales')
plt.plot(Iprueba,recta,'--',color='red',label=r'$\gamma_{c}=\omega_{0}$')
plt.plot(xprueba,recta2,'--',color='red',label=r'$I_{c}$')
plt.xlabel('I(A)')
plt.ylabel(r'$\gamma(s^{-1}$)')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

