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
import ROOT 

from mpl_toolkits.mplot3d import Axes3D


list_delta_L = []
list_random_point = []
list_P_vector = []


###### Dominios de variables #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 10
Theta = np.arange(0, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

##### Tamaño del plano #####
plane_side = 1
long_a = np.arange(-plane_side, plane_side, 0.001)
long_b = np.arange(-plane_side, plane_side, 0.001)

######### Medidas de la CCD ##########
medida_x = 1.197 / 2    # cm
medida_y = 1.587 / 2   # cm
medida_z = 0.0725   # cm

# mapeo_x = np.arange(-medida_x, medida_x, 0.01)
# mapeo_y = np.arange(-medida_y, medida_y, 0.01)
mapeo_z = np.arange(0, medida_z, 0.0001)


mapeo_x = dimension_x()
mapeo_y = dimension_y()
mapeo_z = dimension_z()

# print(mapeo_x)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.


####    Número de Puntos a Simular  ####
number_thet = 10000
number_points_per_angle = 5

n_muons = number_thet * number_points_per_angle

Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

# print(mapeo_x)
#### Inicio del Bucle ####
n_muons_in_CCD = 0

for i in np.arange(0,number_thet):
    list_random_energy_per_plane = []
    list_points_per_plane = []

    Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true

    Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

    Vec = coord_cartesian(Random_th, Random_phi)
    Norma = norma_vec(Vec)
    # print(type(Vec[0]))
    Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
    # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
    # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
    # print('Vector sobre la esfera: ', Point)
    # print('Norma del Vector:', norma)
    # Points.append(Point)

    normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                        ##  hacia el centro de coordenadas
    # normal_Norma_Vec = norma_vec(normal_Vec)
    # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
    # print(len(normal_Vec))
    # Vectors.append(normal_Vec)

    vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
    vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
    # print('Vector Unitario Theta: ', vec_thet)
    # print('Vector Unitario Theta: ', vec_phi)

    for i in np.arange(0,number_points_per_angle):
        flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
        # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
        # list_random_phi.append(Random_phi)


        random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
        random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

        # list_random_a.append(random_a)
        # list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]
        
        list_P_vector.append(P_vector)

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

        # print(random_plane_point)
        list_random_point.append(random_plane_point)

        #### Intersecciones con cada cara   ####
        
        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior ###
        t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (-medida_x - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 4 ###
        t_6 = (medida_x - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
        x_6 = random_plane_point[0] + normal_Vec[0] * t_4

        list_z = [z_3, z_4, z_5, z_6]
        
        if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
            flag_cara_1 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
            flag_cara_2 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
            flag_cara_3 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
            flag_cara_4 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
            flag_cara_5 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
            flag_cara_6 = True
            # print('Bandera 2: ', flag_cara_2)

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        


        Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        if Delta_L != 0:
            if Delta_L > 0:
                list_delta_L.append(Delta_L)
                n_muons_in_CCD = n_muons_in_CCD + muon

            else:
                continue

        else:
                continue
        
Final = datetime.datetime.now()

print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

print('Muones Simulados: ', n_muons)
print('Tamaño de los planos: (' + str(plane_side * 2) + ' x ' + str(plane_side * 2) + ') cm^2')
print('Muones que Impactaron la CCD: ', n_muons_in_CCD)

# fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
# Ax = Axes3D(fig)
# phi = np.arange(0, 2 * np.pi, 0.01)
# theta = np.arange(0, np.pi/2, 0.01)
# theta, phi = np.meshgrid(theta, phi)

# for point in list_random_point:
#     x = point[0]
#     y = point[1]
#     z = point[2]
#     # print(list_points_per_plane[0][i][0])
#     Ax.scatter(x, y, z, c='k', marker='o')

# Ax.scatter(Point[0], Point[1], Point[2], c = 'k', marker='o')

# for point in list_P_vector:
#     x = point[0]
#     y = point[1]
#     z = point[2]
#     # print(list_points_per_plane[0][i][0])
#     Ax.scatter(x, y, z, c='k', marker='o')

# x_s = Radio * np.sin(theta) * np.cos(phi)
# y_s = Radio * np.sin(theta) * np.sin(phi)
# z_s = Radio * np.cos(theta)

# # Agregamos los puntos en el plano 3D
# # ax1.plot_surface(x_s, y_s, z_s)
# Ax.plot_surface(x_s, y_s, z_s)

# Mostramos el gráfico
# plt.show()

array_Delta_L = np.array(list_delta_L)
# max_DeltaL = np.max()
# print()

fig, axs = plt.subplots(figsize=[7,5])
axs.hist(array_Delta_L, bins = 500, label = 'Eventos Simulados: ' + str(number_thet * number_points_per_angle))
# axs.set_xlim(0, 0.5)
axs.vlines([0.0725], 0, 500, colors='k', linestyles='dashed')
# axs.set_xlim(0, 2)
axs.legend()
fig.suptitle(r'Distribución de Longitudes', y = 0.95, size = 20)
plt.show()

file_name = 'array_Delta_L_' + str(n_muons) + '.pkl'
file_object_histogram = open(file_name, 'wb')
pkl.dump(array_Delta_L, file_object_histogram) ## Save the dictionary with all info 
file_object_histogram.close()

h = ROOT.TH1F("h", "Titutlo", 500, 0, 2)
for element in list_delta_L:
    h.Fill(element)


canv = ROOT.TCanvas("canv", "Titulo 2", 700, 400)

canv.cd()
h.Draw()
canv.Print("/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Simulacion_ab_initio/plot.pdf")

# print('Arrray saved in ', current_path + '/' + file_name, ' as a binary file.')
# print('To open use library "pickle". ')