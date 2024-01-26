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

current_path = os.getcwd()

###### Definiciones y dominios #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 100
Theta = np.arange(0, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

long_a = np.arange(0, 10, 0.001)
long_b = np.arange(0, 10, 0.001)

Energy = np.arange(10, 100000, 1)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.

number_thet = 10
number_points_per_angle = 10

n_muons = number_thet * number_points_per_angle

dict_simulation = muon_generator(Energy, Theta, Theta_true, Phi, Radio, long_a, long_b, number_thet, number_points_per_angle)

random_th_array = np.array(dict_simulation['Theta_Radianes'])
random_phi_array = np.array(dict_simulation['Phi_Radianes'])
# type(Random_th)
# len(Random_th)
# Random_th

## Grafica la distribución angular theta ###
fig, axs = plt.subplots(figsize=[7,5])
axs.plot(Theta, 70 * Theta_true)
axs.hist(random_th_array, bins = 110)
fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
plt.show()


## Grafica la distribución angular phi ###
fig, axs = plt.subplots(figsize=[7,5])
axs.hist(random_phi_array, bins = 110)
fig.suptitle(r'Distribución angular $\phi$', y = 0.95, size = 20)
plt.show()

## Gráfica de distribución de los parámetros a y b ##
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=[7,5])
axs[0].hist(dict_simulation['list_random_a'], bins = 110)
axs[1].hist(dict_simulation['list_random_b'], bins = 110)
fig.suptitle(r'Distribución de los parámetros $a$ y $b$', y = 0.95, size = 20)
plt.show()

print('Muones Simulados: ', n_muons)


#### Se crea la figura en 3D
# fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
# Ax = Axes3D(fig)
# plano = 0

# theta = np.arange(0, 2 * np.pi, 0.01)
# phi = np.arange(0, np.pi/2, 0.01)
# theta, phi = np.meshgrid(theta, phi)

# x_s = Radio * np.sin(phi) * np.cos(theta)
# y_s = Radio * np.sin(phi) * np.sin(theta)
# z_s = Radio * np.cos(phi)

# # Agrega la semi-esfera
# # ax1.plot_surface(x_s, y_s, z_s, antialiased=False)
# # Ax.plot_surface(x_s, y_s, z_s, antialiased=False)


# # Toma las coordenadas de cada punto de un plano y las agrega a la figura
# for i in np.arange(0, len(list_points_per_plane[0])):
#     x = list_points_per_plane[plano][i][0]
#     y = list_points_per_plane[plano][i][1]
#     z = list_points_per_plane[plano][i][2]
#     # print(list_points_per_plane[0][i][0])
#     # ax1.scatter(x, y, z, c='k', marker='o')
#     Ax.scatter(x, y, z, c='k', marker='o')

# plt.show()

file_name = 'dict_Simulation_NMUONS_' + str(n_muons) + '.pkl'
file_object_histogram = open(file_name, 'wb')
pkl.dump(dict_simulation, file_object_histogram) ## Save the dictionary with all info 
file_object_histogram.close()

print('Dictionary saved in ', current_path + '/' + file_name, ' as a binary file.')
print('To open use library "pickle". ')
