import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from random_geometry_points.plane import Plane   ### PAra instalar utilizar "pip install random-geometry-points"
from funciones_Sim_ab_initio import *

from mpl_toolkits.mplot3d import Axes3D

###### Definiciones y dominios #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 100
Theta = np.arange(0, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

long_a = np.arange(0, 10, 0.001)
long_b = np.arange(0, 10, 0.001)


Theta_true = dis_angular(Theta) ## Distribución angular theta real.

list_random_th = []
list_random_phi = []
list_random_a = []
list_random_b = []
list_points_per_plane = []
Vectors = []
Points = []

list_random_th = []
list_random_phi = []
list_points_per_plane = []
Vectors = []
Points = []

for i in np.arange(0,len(Theta)):
    list_random_plane_point = []
    Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
    list_random_th.append(Random_th)    ## Lo anexa en una lista

    Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi
    list_random_phi.append(Random_phi)

    # Point = Radio * [ np.sin(Random_th) * np.cos(Random_phi), np.sin(Random_th) * np.sin(Random_phi), np.cos(Random_th)]
    Vec = coord_cartesian(Random_th, Random_phi)
    # print(type(Vec[0]))
    Point = (Radio * Vec[0], Radio * Vec[1], Radio * Vec[2])  ## Genera un punto sobre la esfera.
    Points.append(Point)

    normal_Vec = (-1 * Vec[0], -1 * Vec[1], -1 * Vec[2])     ## Es un vector apuntando hacia el centro de coordenadas
    # print(len(normal_Vec))
    Vectors.append(normal_Vec)

    # plane = Plane(normal_vec = Vec, d_origin = Radio, ref_point = Point, radius = 10.0)   ## Se crea el plano sobre la esfera
    vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
    vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]

    for i in np.arange(0,100):
        # random_plane_point = plane.create_random_points(1)   ## Selección aleatoria de un punto sobre el plano
        # list_random_plane_point.append(random_plane_point)
        random_a = rand.choice(long_a)
        list_random_a.append(random_a)

        random_b = rand.choice(long_b)
        list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], random_a * vec_thet[1] + random_b * vec_phi[1], random_a * vec_thet[2] + random_b * vec_phi[2]]

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        list_random_plane_point.append(random_plane_point)

    list_points_per_plane.append(list_random_plane_point)

random_th_array = np.array(list_random_th)
random_phi_array = np.array(list_random_phi)
# type(Random_th)
# len(Random_th)
# Random_th

## Grafica la distribución angular theta ###
# fig, axs = plt.subplots(figsize=[7,5])
# # # Theta
# axs.plot(Theta, 70 * Theta_true)
# # # plt.plot(Theta_true, Theta)
# # # plt.plot(list_random_th, Theta, '.')
# # # plt.plot(Theta_true, list_random_th, '.')
# axs.hist(random_th_array, bins = 110)
# fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
# plt.show()

## Gráfica de distribución de los parámetros a y b ##
fig, axs = plt.subplots(figsize=[7,5])
axs.hist(list_random_a, bins='auto')
axs.hist(list_random_b, bins = 'auto')
fig.suptitle(r'Distribución angular parámetros $a$ y $b$', y = 0.95, size = 20)
plt.show()


## Phi 
# plt.hist(random_phi_array, bins = 20)
print('Número de planos: ', len(Theta))
print('Número de puntos por plano: ', len(list_points_per_plane[0]))
print('Número total de eventos simulados: ', len(Theta) * len(list_points_per_plane[0]))


#### Se crea la figura en 3D
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
Ax = Axes3D(fig)
plano = 0

theta = np.arange(0, 2 * np.pi, 0.01)
phi = np.arange(0, np.pi/2, 0.01)
theta, phi = np.meshgrid(theta, phi)

x_s = Radio * np.sin(phi) * np.cos(theta)
y_s = Radio * np.sin(phi) * np.sin(theta)
z_s = Radio * np.cos(phi)

# Agrega la semi-esfera
# ax1.plot_surface(x_s, y_s, z_s, antialiased=False)
# Ax.plot_surface(x_s, y_s, z_s, antialiased=False)


# Toma las coordenadas de cada punto de un plano y las agrega a la figura
for i in np.arange(0, len(list_points_per_plane[0])):
    x = list_points_per_plane[plano][i][0]
    y = list_points_per_plane[plano][i][1]
    z = list_points_per_plane[plano][i][2]
    # print(list_points_per_plane[0][i][0])
    # ax1.scatter(x, y, z, c='k', marker='o')
    Ax.scatter(x, y, z, c='k', marker='o')

# Mostramos el gráfico
# plt.show()
