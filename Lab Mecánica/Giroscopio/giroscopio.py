import numpy as np
import matplotlib.pyplot as plt
from math import *
import scipy.optimize as so

from fig_config import (
    add_grid,
    figure_features
)

figure_features()

#Datos del giroscopio y momento de inercia teórico

r = (5.3/2)/100 #r=diam/2
sr = 0.05/100
M= 1.741

I3_teor = M*r**2/2
s_I3_teor = M*r*sr

#Regresion h-t^2 para calcular I3

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

h = np.array([1,0.9,0.8,0.7,0.6])

t_caida_libre = np.array([[5.22,5.19,5.26,5.31,5.29],
        [5.05,5.02,5.01,4.95,5.11],
        [4.65,4.69,4.71,4.70,4.73],
        [4.41,4.42,4.29,4.38,4.28],
        [4.07,4.05,4.18,4.08,4.13]]) 
#array([t100],[t90],[t80],[t70],[t60])

tcaida_media = []
inc_tcaida_media = []

nx,ny = t_caida_libre.shape

for i in range(nx):
    tcaida_media.append(np.mean(t_caida_libre[i]))
    inc_tcaida_media.append(np.std(t_caida_libre[i]))

t_cuadrado = np.array(tcaida_media)**2

reg_h_t = regresionSimple(t_cuadrado,h)

"""fig = plt.figure(figsize=(8,6))
ax = plt.axes(xlim=(15,30),ylim=(0.5,1.05))

ax.set_xlabel(r'$t^2(s^{2})$')
ax.set_ylabel(r'$h(m)$')

add_grid(ax,locations=(0.5,1,1,2))

ax.set_yticks(np.arange(0.5,1.1,0.1))
ax.set_xticks(np.arange(15,31,5))

ax.plot(t_cuadrado,h,'o',color='royalblue',label='Medidas experimentales')
ax.plot(t_cuadrado,reg_h_t[0]+reg_h_t[1]*t_cuadrado,'--',color='navy',label='Recta de regresión')
ax.legend(loc='lower right')
plt.show()"""

b = reg_h_t[1]
sb = reg_h_t[3]

m = 0.06 #masa pequeña en kg
m2 = 0.1
sm = 0.001

I3_reg = (m2*9.8*(r**2))/(2*b)
s_I3_reg = sqrt((sm*(9.8*r**2/(2*b)))**2+(sr*m2*9.8*r/b)**2+(sb*m2*9.8*r**2/(2*b**2))**2)



#Movimiento de precesion

#serie(T,omega)

#Los valores de omega están en medias vueltas, hay que corregirlos

serie1 = np.array([[32.30,28.28,25.28,22.82,20.91],
                [13.7,12.0,10.6,9.6,8.6]])

serie2 = np.array([[38.18,32.51,28.14,24.88,22.66],
                [17.2,14.4,12.6,11.4,10.2]])

serie3 = np.array([[34.72,30.82,27.15,24.59,20.14],
                [15.1,13.1,11.6,9.1,8.3]])

#Para pasar a vueltas por segundo dividimos las frecuencias entre 2 y multiplicamos por 2Pi

serie1[1] = serie1[1]*pi
serie2[1] = serie2[1]*pi
serie3[1] = serie3[1]*pi

s_omega_prec = 0.05
s_T_prec = 0.5

precesion1 = regresionSimple(serie1[1],serie1[0])
precesion2 = regresionSimple(serie2[1],serie2[0])
precesion3 = regresionSimple(serie3[1],serie3[0])

"""fig = plt.figure(figsize=(8,6))
ax = plt.axes(xlim=(25,55),ylim=(15,40))

ax.set_xlabel(r'$\omega(rad\cdot s^{-1})$')
ax.set_ylabel(r'$T(s)$')

add_grid(ax,locations=(0.5,1,1,2))

ax.set_xticks(np.arange(25,55.1,5))
ax.set_yticks(np.arange(15,40.5,5))

ax.plot(serie1[1],precesion1[0]+precesion1[1]*serie1[1],'--',color='b')
ax.plot(serie2[1],precesion2[0]+precesion2[1]*serie2[1],'--',color='g')
ax.plot(serie3[1],precesion3[0]+precesion3[1]*serie3[1],'--',color='r')
ax.plot(serie1[1],serie1[0],marker='.',markersize=10,linestyle='None',color='b',label='Serie 1')
ax.plot(serie2[1],serie2[0],marker='.',markersize=10,linestyle='None',color='g',label='Serie 2')
ax.plot(serie3[1],serie3[0],marker='.',markersize=10,linestyle='None',color='r',label='Serie 3')

ax.legend(loc='lower right')

plt.show()"""

d = 0.067
sd = 0.001

def I_prec(b,sb):
    I3 = (m*9.8*d*b)/(2*pi)
    sI3 = sqrt((sd*m*9.8*b/(2*pi))**2+(sb*m*9.8*d/(2*pi))**2+(sm*9.8*d*b/(2*pi))**2)
    return [I3,sI3]

Iserie1 = I_prec(precesion1[1],precesion1[3])
Iserie2 = I_prec(precesion2[1],precesion2[3])
Iserie3 = I_prec(precesion3[1],precesion3[3])

def media_ponderada(x,sx):
    w = 1/sx**2
    media = np.dot(w,x)/sum(w)
    inc_media = 1/sqrt(sum(w))
    return [media,inc_media]

x = np.array([Iserie1[0],Iserie2[0],Iserie3[0]])
sx = np.array([Iserie1[1],Iserie2[1],Iserie3[1]])

mediaI3 = media_ponderada(x,sx)

print(mediaI3)


###### Obtención de I1

### Muelles

#Cálculo de k

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

    return (sa,sb,sc,media2,excluidos,datos,media1,desviacion1)

st = 0.5

n = 10 #10 oscilaciones

t_k = np.array([4.54,4.41,4.50,4.31,4.3])

inc_t_k = VariasMedidas(t_k)

T = np.array([np.mean(t_k)/10,inc_t_k[2]/10])
T_cuadrado = np.array([T[0]**2,2*T[0]*T[1]])

k = (4*pi**2*m2)/(T[0]**2)
s_k = sqrt(((2*pi**2*m2*T_cuadrado[1])/(T[0]**2))**2+((4*pi**2*sm)/T_cuadrado[0])**2)


#calculo de I1 con muelles

D = 0.325

t_muelles = np.array([8.46,8.61,8.49,8.6,8.52])

inc_t_muelles = VariasMedidas(t_muelles)
T_muelles = np.array([np.mean(t_muelles)/10,inc_t_muelles[2]/10]) #T,sT

omega = np.array([2*pi/T_muelles[0],2*pi*T_muelles[1]/T_muelles[0]**2])

I1_muelle = 2*k*(D**2)/(omega[0]**2)
s_I1_muelle = sqrt((2*s_k*D/omega[0]**2)**2+(4*D*k*0.001/omega[0]**2)**2+(k*D*omega[1]/omega[0]**3)**2)

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



"""fig = plt.figure(figsize=(8,6))
ax = plt.axes(xlim=(15,55),ylim=(0,9))

ax.set_xlabel(r'$\omega(rad\cdot s^{-1})$')
ax.set_ylabel(r'$\Omega_{n}(rad\cdot s^{-1})$')

add_grid(ax,locations=(0.5,1,1,2))

ax.set_xticks(np.arange(15,55.1,5))
ax.set_yticks(np.arange(0,9.1,1))

ax.errorbar(omega1[0],omega_nut1[0],yerr=omega_nut1[1],linestyle='None',marker='.',markersize=10,color='b',label='Serie 1')
ax.plot(omega1[0],reg_nut1[0]+reg_nut1[1]*omega1[0],'--',color='b')
ax.errorbar(omega2[0],omega_nut2[0],yerr=omega_nut2[1],linestyle='None',marker='.',markersize=10,color='g',label='Serie 2')
ax.plot(omega2[0],reg_nut2[0]+reg_nut2[1]*omega2[0],'--',color='g')
ax.errorbar(omega3[0],omega_nut3[0],yerr=omega_nut3[1],linestyle='None',marker='.',markersize=10,color='r',label='Serie 3')
ax.plot(omega3[0],reg_nut3[0]+reg_nut3[1]*omega3[0],'--',color='r')

ax.legend(loc='lower right')

plt.show()"""

#### I1 teorico

m_disco = 1.741
m_contrapeso = 0.895

sm = 0.0001
s_d = 0.001

d_disco = 10.5/100
ds = np.array([20.5,23.5])/100
d_contrapeso = np.mean(ds)
s_d_contrapeso = np.std(ds)

I1_teor = m_disco*d_disco**2 + m_contrapeso*d_contrapeso**2
s_I1_teor = sqrt((sm*d_disco**2)**2+(sm*d_contrapeso**2)**2+(2*m_disco*d_disco*s_d)**2+(2*m_contrapeso*d_contrapeso*s_d_contrapeso)**2)







