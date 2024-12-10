import numpy as np
from math import *
from numpy.linalg import *

import scipy.optimize as so

import matplotlib.pyplot as plt

from sympy.solvers import solve
from sympy import Symbol

def VariasMedidas(x):   
    k = 2 #Factor de cobertura
    
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
    
    return (media2,sc)

### I = 0,3A

sb = 0.5 #Amplitud de tipo b para el tiempo

v2 = np.array([21.86,21.96,22.08])
v4 = np.array([13.10,13.52,13.17,13.17])
v6 = np.array([8.97,8.84,8.86,8.9])
v8 = np.array([6.54,6.35,6.63,6.56])
v10 = np.array([5.39,5.31,5.25,5.38])
v12 = np.array([4.54,4.41,4.40,4.41])
v9 = np.array([5.91,5.91,5.90,5.86])
v7 = np.array([7.48,7.88,7.73,7.54])
v6_5 = np.array([8.26,8.2,8.2,8.26])
v7_5 = np.array([7.15,7.15,7.08,7.08])
v7_25 = np.array([7.4,7.47,7.55,7.4])
v6_75 = np.array([8.19,8.11,8.07,8.13])

def media(t):
    array = VariasMedidas(t)
    print(f'Media: {array[0]} y su inc: {array[1]}')
    return array

def frecuencia(t,st,n):
    omega1 = 2*pi*n/t
    T = 2*pi/omega1
    s_omega1 = st/(t**2)
    print('Para t=',t,' la frecuencia es ',omega1,' y la incertidumbre es ',s_omega1,'El periodo es: ',T)
    return(omega1,s_omega1,T)

sb=0.2 #Amplitud de tipo b para la amplitud

t2=media(v2)
omega2=frecuencia(t2[0],t2[1],4)
t4=media(v4)
omega4=frecuencia(t4[0],t4[1],4)
t6=media(v6)
omega6=frecuencia(t6[0],t6[1],4)
t8=media(v8)
omega8=frecuencia(t8[0],t8[1],4)
t10=media(v10)
omega10=frecuencia(t10[0],t10[1],4)
t12=media(v12)
omega12=frecuencia(t12[0],t12[1],4)
t9=media(v9)
omega9=frecuencia(t9[0],t9[1],4)
t7=media(v7)
omega7=frecuencia(t7[0],t7[1],4)
t6_5=media(v6_5)
omega6_5=frecuencia(t6_5[0],t6_5[1],4)
t7_5=media(v7_5)
omega7_5=frecuencia(t7_5[0],t7_5[1],4)
t7_25=media(v7_25)
omega7_25=frecuencia(t7_25[0],t7_25[1],4)
t6_75=media(v6_75)
omega6_75=frecuencia(t6_75[0],t6_75[1],4)

omega_fuerza3 = np.array([omega2[0],omega4[0],omega6[0],omega8[0],omega10[0],omega12[0],omega9[0],omega7[0],omega6_5[0],omega7_5[0],omega7_25[0],omega6_75[0]])

a2 = [0.4,0.6]
a4 = [0.6,0.6]
a6 = [1.4,1.2]
a8 = [1.6,1.4]
a10 = [0.4,0.6]
a12 = [0.2,0.4]
a9 = [0.8,0.6]
a7 = [7.4,7.0]
a6_5 = [2.2,2.4]
a7_5 = [2.8,2.6]
a7_25 = [5.2,4.8]
a6_75 = [3.2,2.8]

am2=media(a2)
am4=media(a4)
am6=media(a6)
am8=media(a8)
am10=media(a10)
am12=media(a12)
am9=media(a9)
am7=media(a7)
am6_5=media(a6_5)
am7_5=media(a7_5)
am7_25=media(a7_25)
am6_75=media(a6_75)

amplitudesI3 = np.array([am2[0],am4[0],am6[0],am8[0],am10[0],am12[0],am9[0],am7[0],am6_5[0],am7_5[0],am7_25[0],am6_75[0]])

s_amplitudesI3 = np.array([am2[1],am4[1],am6[1],am8[1],am10[1],am12[1],am9[1],am7[1],am6_5[1],am7_5[1],am7_25[1],am6_75[1]])



##### I=0,6 A

sb = 0.5 #Amplitud de tipo b para el tiempo

V2 = np.array([32.65,32.52,32.46,32.42])
V4 = np.array([13.62,13.57,13.56,13.58])
V6 = np.array([8.98,8.99,8.92,8.91])
V8 = np.array([6.64,6.7,6.69,6.76])
V10 = np.array([5.26,5.32,5.31,5.45])
V12 = np.array([4.34,4.46,4.39,4.47])
V7 = np.array([7.55,7.6,7.67,7.67])
V6_5 = np.array([8.21,8.13,8.18,8.26])
V7_5 = np.array([7.28,7.21,7.22,7.2])
V7_25 = np.array([7.34,7.41,7.54,7.48])

t2=media(V2)
omega2=frecuencia(t2[0],t2[1],4)
t4=media(V4)
omega4=frecuencia(t4[0],t4[1],4)
t6=media(V6)
omega6=frecuencia(t6[0],t6[1],4)
t8=media(V8)
omega8=frecuencia(t8[0],t8[1],4)
t10=media(V10)
omega10=frecuencia(t10[0],t10[1],4)
t12=media(V12)
omega12=frecuencia(t12[0],t12[1],4)
t7=media(V7)
omega7=frecuencia(t7[0],t7[1],4)
t6_5=media(V6_5)
omega6_5=frecuencia(t6_5[0],t6_5[1],4)
t7_5=media(V7_5)
omega7_5=frecuencia(t7_5[0],t7_5[1],4)
t7_25=media(V7_25)
omega7_25=frecuencia(t7_25[0],t7_25[1],4)

omega_fuerza6 = np.array([omega2[0],omega4[0],omega6[0],omega8[0],omega10[0],omega12[0],omega7[0],omega6_5[0],omega7_5[0],omega7_25[0]])

sb = 0.2 #Inc de la amplitud

A2 = [0.4,0.6]
A4 = [0.6,0.6]
A6 = [1.2,1]
A8 = [1.2,1.4]
A10 = [0.4,0.6]
A12 = [0.2,0.4]
A7 = [2,1.8]
A6_5 = [1.6,1.4]
A7_5 = [1.8,1.6]
A7_25 = [2,1.8]


Am2=media(A2)
Am4=media(A4)
Am6=media(A6)
Am8=media(A8)
Am10=media(A10)
Am12=media(A12)
Am7=media(A7)
Am6_5=media(A6_5)
Am7_5=media(A7_5)
Am7_25=media(A7_25)

amplitudesI6 = np.array([Am2[0],Am4[0],Am6[0],Am8[0],Am10[0],Am12[0],Am7[0],Am6_5[0],Am7_5[0],Am7_25[0]])

s_amplitudesI6 = np.array([Am2[1],Am4[1],Am6[1],Am8[1],Am10[1],Am12[1],Am7[1],Am6_5[1],Am7_5[1],Am7_25[1]])

def amplitud(omega,F0_J,omega0,gamma):
    y = F0_J/(((omega**2-omega0**2)**2+4*(gamma**2)*(omega**2))**(1/2))
    return y

def ajuste(fun,par,x,y):
    sol = so.curve_fit(fun,x,y,p0=(par)) #Array de soluciones y incertidumbres
    valores = sol[0]
    inc = np.sqrt(np.diag(sol[1]))
    print(f'Los valores de los parametros a ajustar son {valores}')
    print(f'Las incertidumbres son: {inc}')
    return valores

par = [0.1,0.1,0.1]

amplitud_w_3 = ajuste(amplitud,par,omega_fuerza3,amplitudesI3)
amplitud_w_6 = ajuste(amplitud,par,omega_fuerza6,amplitudesI6)


xprueba = np.arange(0,8,0.0001)
yEst1 = amplitud(xprueba,amplitud_w_3[0],amplitud_w_3[1],amplitud_w_3[2])
yEst2 = amplitud(xprueba,amplitud_w_6[0],amplitud_w_6[1],amplitud_w_6[2])

a_r3 = np.where(yEst1 == np.amax(yEst1))
a_r6 = np.where(yEst2 == np.amax(yEst2))

omega_r3 = xprueba[a_r3]
omega_r6 = xprueba[a_r6]

gamma3 = 0.10221
s_gamma3 = 0.00092
gamma6 = 0.353
s_gamma6 = 0.011
omega_0 = 3.32958
s_omega_0 = 0.00016

inc_comparacion3 = sqrt((2*omega_0*s_omega_0)**2+(4*gamma3*s_gamma3)**2)
inc_comparacion6 = sqrt((2*omega_0*s_omega_0)**2+(4*gamma6*s_gamma6)**2)

print(f'omega_r^2 para 0,3: {omega_r3**2} y para 0,6: {omega_r6**2}')
print(f'Los valores de omega_0^2-gamma^2 son {omega_0**2-2*gamma3**2} y {omega_0**2-2*gamma6**2} y sus incertidumbres {inc_comparacion3} y {inc_comparacion6}')



"""plt.errorbar(omega_fuerza3,amplitudesI3,s_amplitudesI3,marker='o',color='black',label='I=0,3 A',linestyle='None')
plt.errorbar(omega_fuerza6,amplitudesI6,s_amplitudesI6,marker='o',color='grey',label='I=0,6 A',linestyle='None')
plt.plot(xprueba,yEst1)
plt.plot(xprueba,yEst2,color='red')
plt.xlabel(r'$\omega(s^{-1})$')
plt.ylabel(r'$\theta_{max}$')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()"""

#Factor de calidad

amax_3 = np.amax(amplitudesI3)
amax_6 = np.amax(amplitudesI6)

altura_anchura_3 = amax_3/sqrt(2)
altura_anchura_6 = amax_6/sqrt(2)

anchura3 = 2*gamma3
anchura6 = 2*gamma6

omega_r = omega_0 #Para gamma muy pequeñaas

def factorQ(omega_r,anchura,s_omega_r,s_gamma):
    Q = omega_r/anchura
    sQ = sqrt((s_omega_r/(anchura))**2+((s_gamma*omega_r)/(2*(anchura/2)**2))**2)
    return [Q,sQ]

Q3 = factorQ(omega_r,anchura3,s_gamma3,s_omega_0)
Q6 = factorQ(omega_r,anchura6,s_gamma6,s_omega_0)

print(f'Para I=0,3 Q: {Q3} y para 0,6 Q: {Q6}')

#### Grafica del desfase

def desfase(x,gamma,omega0):
    y = np.arctan(((2*gamma*x)/((omega0**2)-(x**2))))
    return y


desfase3 = desfase(xprueba,gamma3,omega_0)
for i in range(len(desfase3)):
    if desfase3[i] < 0:
        desfase3[i] += pi

desfase6 = desfase(xprueba,gamma6,omega_0)
for i in range(len(desfase6)):
    if desfase6[i] < 0:
        desfase6[i] += pi



"""plt.plot(xprueba,desfase3,color='royalblue',label='I=0,3 A')
plt.plot(xprueba,desfase6,color='red',label='I=0,6 A')
plt.xlabel(r'$\omega (s^{-1})$')
plt.ylabel(r'$\delta(rad)$')
plt.grid(True)
plt.legend(loc='upper left')
plt.show()"""