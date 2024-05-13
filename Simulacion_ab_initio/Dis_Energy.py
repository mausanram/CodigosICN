import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
# from random_geometry_points.plane import Plane   ### Para instalar utilizar "pip install random-geometry-points"
from funciones_Sim_ab_initio import *
import datetime
import os
import pickle as pkl
import scipy.integrate as integrate

### Configuración de estilo de las gráficas ###
plt.rcParams.update({
    "image.origin": "lower",
    "image.aspect": 1,
    #"text.usetex": True,
    "grid.alpha": .5,
    "axes.linewidth":2,
    "lines.linewidth" : 1,
    "font.size":    15.0,
    "xaxis.labellocation": 'right',  # alignment of the xaxis label: {left, right, center}
    "yaxis.labellocation": 'top',  # alignment of the yaxis label: {bottom, top, center}
    "xtick.top":           True ,  # draw ticks on the top side
    "xtick.major.size":    8    ,# major tick size in points
    "xtick.minor.size":    4      ,# minor tick size in points
    "xtick.direction":     'in',
    "xtick.minor.visible": True,
    "ytick.right":           True ,  # draw ticks on the top side
    "ytick.major.size":    8    ,# major tick size in points
    "ytick.minor.size":    4      ,# minor tick size in points
    "ytick.direction":     'in',
    "ytick.minor.visible": True,
    "ytick.major.width":   2   , # major tick width in points
    "ytick.minor.width":   1 ,
    "xtick.major.width":   2   , # major tick width in points
    "xtick.minor.width":   1 ,
    "legend.framealpha": 0 ,
    "legend.loc": 'best',
})


Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

list_Energy = []
# E = np.arange(10, 100000, 100)
E_in = 10 ** (-1)
E_fin = 10 ** 6
N = 10000

list_Energy = Energy_list(E_in, E_fin, N)
print('Longitud de la lista de energía: ', len(list_Energy))


# Thet = [0.00174533, 0.523599, 0.785398, 1.309]     ## En radianes
# Ang = [0, 30, 45, 75]   ## Grados

Thet = [0.00174533, 0.785398, 1.309, np.radians(80)]
Ang = [0, 45, 75, 80]

# Thet = [0.0174533, 1.0472,  0.785398, 1.309 , 1.55334]   
# Ang = [0, 45, 60, 75, 90]

# Thet = [0.0000000000000001] 
# Ang = [0]

# for element in np.arange(0, len(Thet)):
# # for element in np.arange(0, 1):
#     Energy = dis_energy(E, Thet[element])
#     # print(Energy)
#     plt.plot(E, Energy, label = str(Ang[element]) + '°')
fig, axs = plt.subplots(1,2,figsize=[15,5])

for element in np.arange(0, len(Thet)):
    list_dis_Energy = []
    for energy in list_Energy:
        Energy = dis_energy(energy, Thet[element])
        list_dis_Energy.append(Energy)
        # print(Energy)
    axs[0].plot(list_Energy, list_dis_Energy, label = r'$\theta = $' + str(Ang[element]) + '°')

Ang = [0.01, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 85, 90]
Ang_array = np.arange(0, 90, 1)
Thet = []
Thet_array = []
list_results = []

for ang in Ang:
    rad = np.radians(ang)
    Thet.append(rad)

for ang in Ang_array:
    rad = np.radians(ang)
    Thet_array.append(rad)

for angle in Thet: 
    result = integrate.quad(dis_energy, a = 0, b = np.inf, args = angle)
    list_results.append(result[0]/ (3.52877403746463 *  10**(-5) ))
    del result

res = integrate.quad(dis_energy, a = 0, b = np.inf, args = np.radians(80))
Res = res[0] / (3.52877403746463 *  10**(-5) )
Real = np.cos(np.radians(80))**2

axs[1].plot(Ang, list_results, 'ob' )
axs[1].plot(80, Real - Res, 'ob')
axs[1].plot(Ang_array, np.cos(Thet_array)**2, 'k')
axs[1].set_xlabel('Ángulo (°)')
axs[1].set_ylabel('I / I_0')
axs[1].grid()


axs[0].set_xlabel('Energy (MeV)')   
axs[0].grid() 
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend()
axs[0].set_title('Distribuciónes de Smith-Duller')

plt.show()

Final = datetime.datetime.now()

print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

del Energy 