# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 19:20:26 2023

@author: User
"""

from numpy import *
from numpy.linalg import *
from pylab import *
from numpy.polynomial.polynomial import Polynomial
import matplotlib.pyplot
def pondreg(x,y,s):

    w=1.0/(s*s)
    wy=sum(w*y); wx=sum(w*x);
    wxx=sum(w*x*x); wxy=sum(w*x*y); wyy=sum(w*y*y)
    sw=sum(w)
    d=det(array([[sw, wx],[wx, wxx]]))
    a=(wy*wxx-wx*wxy)/d
    b=(sw*wxy-wx*wy)/d
    sa=sqrt(wxx/d); sb=sqrt(sw/d)
    r=(sw*wxy-wx*wy)/sqrt((sw*wxx-wx**2)*(sw*wyy-wy**2))
    
    p1=a+b*min(x)
    p2=a+b*max(x)

    xp=array((min(x),max(x)))
    yp=array((p1,p2))
    
    return plot(xp,yp,'r-',label='Ajuste')
def simpreg(x,y):

    n=len(x)
    sx=sum(x); sy=sum(y); xx=dot(x,x); yy=dot(y,y); xy=dot(x,y);
    denom=(n*xx - sx**2)
    b=(n*xy - sx*sy)/denom
    a=(xx*sy - sx*xy)/denom
    s=sqrt(sum((y-a-b*x)**2)/(n-2))
    sa=s*sqrt(xx/(n*xx-sx**2))
    sb=s*sqrt(n/(n*xx-sx**2))
    r=(n*xy-sx*sy)/sqrt((n*xx-sx**2)*(n*yy-sy**2))
    
    p1=a+b*min(x)
    p2=a+b*max(x)

    xp=array((min(x),max(x)))
    yp=array((p1,p2))
    
    return plot(xp,yp,'r-',label='Ajuste')

x = array((0.118952724, 0.203076124, 0.28517594, 0.369750319, 0.456251728,0.544498912, 0.633058306, 0.723500761, 0.81403771, 0.906402816,1))
v = array((14.52
,15.68
,15.99
,16.12
,16.22
,16.35
,16.55
,16.83
,17.18
,17.609
,18.117))
s=array((0.65
,0.39
,0.28
,0.23
,0.19
,0.17
,0.15
,0.13
,0.11
,0.095
,0.080))
figure(dpi=200)
pondreg(x,v,s)
plot(x,v,'b.')
errorbar(x,v,yerr=s, capsize=2, linestyle='None',color='blue')
grid('True')
xlabel('$\chi_{H_{2}O}$')
ylabel('$v_{H_{2}O}(cm^{3}/mol)$')

show()