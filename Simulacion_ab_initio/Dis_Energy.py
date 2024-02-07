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

E = np.arange(10, 100000, 100)
Thet = [0.00174533, 0.523599, 0.785398, 1.309]     ## En radianes
Ang = [0, 30, 45, 75]

for element in np.arange(0, len(Thet)):
# for element in np.arange(0, 1):
    Energy = dis_energy(E, Thet[element])
    # print(Energy)
    plt.plot(E, Energy, label = str(Ang[element]) + '°')

plt.xlabel('Energy (MeV)')   
plt.grid() 
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

Final = datetime.datetime.now()

print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

del Energy 
del E