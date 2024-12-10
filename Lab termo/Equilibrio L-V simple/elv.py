import numpy as np
import matplotlib.pyplot as plt

lnp = np.array([8.517193191, 8.5754621, 8.612503371, 8.809862805, 8.853665428, 8.9226583, 9.024010794, 9.148464968, 9.26860928, 9.384293679, 9.503009986, 9.642122788, 9.752664663, 9.857443614, 9.98506754, 10.07744086, 10.19988157, 10.30226367, 10.39818379, 10.49127422, 10.58910647, 10.68051622, 10.76214917, 10.84544603, 10.94022444, 11.02841715, 11.10945836, 11.19409666, 11.27466828, 11.35392973, 11.44035477, 11.51990108, 11.60093634, 11.67674355, 11.75194237, 11.83500896, 11.91438255, 11.99226042, 12.06220228, 12.13618652, 12.21106019, 12.28025204, 12.34930681, 12.41633357, 12.4848401, 12.55602504, 12.62445398, 12.68788132, 12.74606865, 12.81338957, 12.8741584, 12.94203982, 12.99634654, 13.05814481, 13.12016095])


T = np.array([0.003448276, 0.003424658, 0.003401361, 0.003378378, 0.003355705, 0.003333333, 0.003311258, 0.003289474, 0.003267974, 0.003246753, 0.003225806, 0.003205128, 0.003184713, 0.003164557, 0.003144654, 0.003125, 0.00310559, 0.00308642, 0.003067485, 0.00304878, 0.003030303, 0.003012048, 0.002994012, 0.00297619, 0.00295858, 0.002941176, 0.002923977, 0.002906977, 0.002890173, 0.002873563, 0.002857143, 0.002840909, 0.002824859, 0.002808989, 0.002793296, 0.002777778, 0.002762431, 0.002747253, 0.00273224, 0.002717391, 0.002702703, 0.002688172, 0.002673797, 0.002659574, 0.002645503, 0.002631579, 0.002617801, 0.002604167, 0.002590674, 0.00257732, 0.002564103, 0.00255102, 0.002538071, 0.002525253, 0.002512563])

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
    sx=np.sum(x); sy=np.sum(y); xx=np.dot(x,x); yy=np.dot(y,y); xy=np.dot(x,y);
    denom=(n*xx - sx**2)
    b=(n*xy - sx*sy)/denom
    a=(xx*sy - sx*xy)/denom
    s=np.sqrt(sum((y-a-b*x)**2)/(n-2))
    sa=s*np.sqrt(xx/(n*xx-sx**2))
    sb=s*np.sqrt(n/(n*xx-sx**2))
    r=(n*xy-sx*sy)/np.sqrt((n*xx-sx**2)*(n*yy-sy**2))
    return [a,b, sa, sb, r, s]



reg = regresionSimple(T,lnp)

l = -reg[1]*8.31
print(l)























