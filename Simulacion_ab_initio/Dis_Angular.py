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

from mpl_toolkits.mplot3d import Axes3D

###### Definiciones y dominios #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 100
Theta = np.arange(0, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

long_a = np.arange(-5, 5, 0.01)
long_b = np.arange(-5, 5, 0.01)

# Energy = np.arange(10, 100000, 1)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.

number_thet = 60000
number_points_per_angle = 1

n_muons = number_thet * number_points_per_angle

list_random_thet = []
list_random_phi = []

for i in np.arange(0,number_thet):
    list_points_per_plane = []
    list_random_energy_per_plane = []

    Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
    list_random_thet.append(Random_th)

    Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi
    list_random_phi.append(Random_phi)

    # Vec = coord_cartesian(Random_th, Random_phi)
    # Point = (Radio * Vec[0], Radio * Vec[1], Radio * Vec[2])  ## Genera un punto sobre la esfera.

Bins = 110
delta_B = 1.6 / Bins
Const_Normal =  3 * number_thet * delta_B
print('La constante de ajuste de la curva es: ', np.round(Const_Normal, 3))


fig, axs = plt.subplots(figsize=[7,5])
axs.plot(Theta, Const_Normal * Theta_true)
axs.hist(np.array(list_random_thet), bins = Bins)
fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
plt.show()


fig, axs = plt.subplots(figsize=[7,5])
# axs.plot(Theta, Const_Normal * Theta_true)
axs.hist(list_random_phi, bins = Bins)
fig.suptitle(r'Distribución angular $\phi$', y = 0.95, size = 20)
plt.show()
