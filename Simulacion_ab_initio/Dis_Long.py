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


list_delta_L = []
list_random_point = []
list_P_vector = []


###### Definiciones y dominios #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 100
Theta = np.arange(0, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

long_a = np.arange(-5, 5, 0.01)
long_b = np.arange(-5, 5, 0.01)

# Energy = np.arange(10, 100000, 1)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.

number_thet = 20000
number_points_per_angle = 1

n_muons = number_thet * number_points_per_angle

for i in np.arange(0,number_thet):
        list_points_per_plane = []
        list_random_energy_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true

        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        # Energy_dist = dis_energy(E, Random_th[0])

        Vec = coord_cartesian(Random_th, Random_phi)
        # print(type(Vec[0]))
        Point = (Radio * Vec[0], Radio * Vec[1], Radio * Vec[2])  ## Genera un punto sobre la esfera.
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0], -1 * Vec[1], -1 * Vec[2])     ## Es un vector apuntando hacia el centro de coordenadas
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]

        for i in np.arange(0,number_points_per_angle):
            # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            # list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # list_random_a.append(random_a)
            # list_random_b.append(random_b)

            delta_L = 0.0725 / np.cos(Random_th)  ## cm
            list_delta_L.append(delta_L)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            list_P_vector.append(P_vector)

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            list_random_point.append(random_plane_point)

            # Random_energy = rand.choices(E, Energy_dist)
            # list_random_energy.append(Random_energy[0])

# fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
# Ax = Axes3D(fig)
# theta = np.arange(0, 2 * np.pi, 0.01)
# phi = np.arange(0, np.pi/2, 0.01)
# theta, phi = np.meshgrid(theta, phi)

# for point in list_random_point:
#     x = point[0]
#     y = point[1]
#     z = point[2]
#     # print(list_points_per_plane[0][i][0])
#     Ax.scatter(x, y, z, c='k', marker='o')

# for point in list_P_vector:
#     x = point[0]
#     y = point[1]
#     z = point[2]
#     # print(list_points_per_plane[0][i][0])
#     Ax.scatter(x, y, z, c='k', marker='o')

# x_s = Radio * np.sin(phi) * np.cos(theta)
# y_s = Radio * np.sin(phi) * np.sin(theta)
# z_s = Radio * np.cos(phi)

# Agregamos los puntos en el plano 3D

# ax1.plot_surface(x_s, y_s, z_s)

# Mostramos el gráfico
plt.show()

fig, axs = plt.subplots(figsize=[7,5])
axs.hist(np.array(list_delta_L), bins = 400)
axs.set_xlim(0, 0.5)
axs.vlines([0.0725], 0, 600, colors='k', linestyles='dashed')
fig.suptitle(r'Distribución del Delta L', y = 0.95, size = 20)
plt.show()