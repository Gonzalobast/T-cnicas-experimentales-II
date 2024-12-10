import numpy as np
from math import *
import matplotlib.pyplot as plt
from numpy.linalg import *
import scipy.optimize as so
import pandas as pd


#Disposicion 1

t1 = []
for i in range(238):
    t1.append(10*i)

datos1 = pd.read_csv('Balanza de torsión\datos1.csv',sep=';')
datos2 = pd.read_csv('Balanza de torsión\datos2.csv',sep=';')

t1 = np.array(datos1.iloc[:,0])
x1 = np.array(datos1.iloc[:,1])
t2 = np.array(datos2.iloc[:,0])
x2 = np.array(datos2.iloc[:,1])

sx = 1
st = 0.5

def x(t,x0,A,b,c,d):
    x =x0+A*np.exp(-b*t)*np.cos(c*t+d)
    return x

par = [20,50,0.01,0.01,0.1]

sol_1 = so.curve_fit(x,t1,x1,p0=(par))
sol_2 = so.curve_fit(x,t2,x2,p0=(par))

#Valores de los coeficientes y incertidumbres
coeff1 = sol_1[0]
inc1 = np.diag(sol_1[1])

coeff2 = sol_2[0]
inc2 = np.diag(sol_1[1])

curva1 = x(t1,coeff1[0],coeff1[1],coeff1[2],coeff1[3],coeff1[4])
curva2 = x(t2,coeff2[0],coeff2[1],coeff2[2],coeff2[3],coeff2[4])

x0_1 = []
x0_2 = []

for i in range(len(x1)):
    x0_1.append(coeff1[0])

for i in range(len(x2)):
    x0_2.append(coeff2[0])

omega1 = coeff1[3]
s_omega1 = inc1[3]
omega2 = coeff2[3]
s_omega2 = inc2[3]

def periodo(omega,s_omega):
    T = 2*pi/omega
    sT = (2*pi*s_omega)/(omega**2)
    return [T,sT]

periodo1 = periodo(omega1,s_omega1)
periodo2 = periodo(omega2,s_omega2)

w_a = 1/(periodo1[1]**2)
w_b = 1/(periodo2[1]**2)

T = (w_a*periodo1[0]+w_b*periodo2[0])/(w_a+w_b) #Media ponderada
incT = 1/(sqrt(w_a+w_b))

print(f'El periodo es: {T} y la incertidumbre es: {incT}')

plt.plot(t2,x2,'.',color='orange',label='Medidas experimentales')
plt.plot(t2,curva2,'--',color='black',label='Ajuste')
plt.plot(t2,x0_2,'--',color='red',label=r'$x_{eq}$')
plt.xlabel('t(s)')
plt.ylabel('x(cm)')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()

#Calculo de los residuos (Diferencia entre la curva y los valores medidos)

residuos1 = abs(curva1-x1)
residuos2 = abs(curva2-x2)

"""plt.errorbar(t2,residuos2,sx,color='orange',marker='.',linestyle='None',ecolor='grey',label=r'$x_{err}$ y su incertidumbre')
plt.xlabel('t(s)')
plt.ylabel(r'$x_{err}(cm)$')
plt.grid(True)
plt.legend(loc='upper left')
plt.show()"""

#Calculo de la distancia L

l1 = 15
s_l1 = 0.1
l2 = 22
s_l2 = s_l1
l = 560
sl=1.4
l_a=17.1
s_la=1.5

L = (l1+l2+l-l_a)/100
s_L = (s_l1*2+sl+s_la)/100 #pasamos a metros

#calculo de b

d = 0.05 #distancia entre el eje de giro y las bolas en m 

"""def b(xeq):
    b = ((xeq)*d)/(2*(L))
    sb = sqrt(((sx/100)*d/(2*(L)))**2+((s_L)*(xeq)*d/(2*(L)**2))**2)
    return [b,sb]

b1 = b(coeff1[0]/100) #hay que pasar a metros
b2 = b(coeff2[0]/100)

w_b1 = 1/(b1[1]**2)
w_b2 = 1/(b2[1]**2)

b = (w_b1*b1[0]+w_b2*b2[0])/(w_b1+w_b2)
s_b = 1/sqrt(w_b1+w_b2)"""

M=1.5
sM = 0.01
delta_x = abs(coeff1[0]-coeff2[0])/100
s_deltax = 2*sx/100
Gn = 6.67e-11

#s-deltax

s_deltax = (inc1[0]+inc2[0])/100

b = sqrt((6.67e-11*M*T**2*L)/(delta_x*pi**2*d))
s_bcuadrado = sqrt(((sM*Gn*T**2*L)/(delta_x*pi**2*d))**2+((2*T*incT*Gn*M*L)/(delta_x*pi**2*d))**2+((s_L*Gn*T**2*M)/(delta_x*pi**2*d))**2+((s_deltax*Gn*M*T**2*L)/(delta_x**2*pi**2*d))**2)
s_b = s_bcuadrado/(2*b)

print(f'b= {b} sb= {s_b}')

#Calculo de G

beta = b**3/((b**2+4*d**2)**(3/2))
s_beta = s_b*((3*b**2)/(4*d**2+b**2)**(3/2)-(3*b**4)/(4*d**2+b**2)**(5/2))

G = (pi**2*b**2*d*delta_x)/(M*(L)*T**2*(1-beta))

print(G)

sG = sqrt(((2*b*s_b*pi**2*d*delta_x)/(M*L*T**2*(1-beta)))**2+((s_deltax*pi**2*b**2*d)/(M*L*T**2*(1-beta)))**2+((sM*pi**2*b**2*d*delta_x)/(M**2*L*T**2*(1-beta)))**2+((s_L*pi**2*b**2*d*delta_x)/(M*L**2*T**2*(1-beta)))**2+((incT*pi**2*b**2*d*delta_x*2)/(M*L*T**3*(1-beta)))**2+((s_beta*pi**2*b**2*d*delta_x)/(M*L*T**2*(1-beta)**2))**2)

print(f'sG: {sG}')

#Analisis de las incertidumbres y los aportes de cada termino

ub = (2*b*s_b*pi**2*d*delta_x)/(M*L*T**2*(1-beta))
udelta_x = (s_deltax*pi**2*b**2*d)/(M*L*T**2*(1-beta))
uM = (sM*pi**2*b**2*d*delta_x)/(M**2*L*T**2*(1-beta))
uL = (s_L*pi**2*b**2*d*delta_x)/(M*L**2*T**2*(1-beta))
uT = (incT*pi**2*b**2*d*delta_x*2)/(M*L*T**3*(1-beta))
ubeta = (s_beta*pi**2*b**2*d*delta_x)/(M*L*T**2*(1-beta)**2)

print(f'ub: {ub}, udeltax: {udelta_x}, uM: {uM}, uL: {uL}, uT:{uT}, ubeta:{ubeta}')

# Calculo de la constante de torsion

rM = 0.03158 #Radios en metros
rm = 0.0082

dens = M/((4/3)*pi*(rM**3))
inc_dens = sM/(4/3*pi*rM**3)

m = dens*(4/3*pi*(rm**3))
sm = inc_dens*(4/3*pi*rm**3)

k = (8*pi**2*m*d**2)/T**2

s_k = sqrt((sm*pi**2*8*d**2/T**2)**2+(incT*16*pi**2*m*d**2/T**3)**2)

print(f'k:{k} y inc:{s_k}')