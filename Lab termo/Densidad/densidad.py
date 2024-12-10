import matplotlib.pyplot as plt
import numpy as np
from math import * 

#Datos de la memoria

rho_96 = 0.805
rho_agua = 1
rho_oh = 0.789

m_96 = np.array([2.019, 3.939, 5.85, 7.749, 9.724, 11.752, 13.802, 15.826, 17.806, 20.403, 0])

m_agua = np.array([18.001, 16.101, 14.15, 12.314, 10.749, 8.408, 6.241, 4.238, 2.206, 0, 20.211])

m_total = m_96 + m_agua

m_oh = 0.96*m_96*rho_oh/rho_96
m_agua_total = m_agua + 0.04*m_96/rho_96

