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

Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

list_Energy = []
# E = np.arange(10, 100000, 100)
E_in = 10 ** (-1)
E_fin = 10 ** 6
N = 1000

list_Energy = Energy_list(E_in, E_fin, N)


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

for element in np.arange(0, len(Thet)):
    list_dis_Energy = []
    for energy in list_Energy:
        Energy = dis_energy(energy, Thet[element])
        list_dis_Energy.append(Energy)

        # print(Energy)
    plt.plot(list_Energy, list_dis_Energy, label = str(Ang[element]) + '°')

plt.xlabel('Energy (MeV)')   
plt.grid() 
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.title('Distribución de Energía (Smith-Duller)')
plt.show()

Final = datetime.datetime.now()

print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

del Energy 