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
# import ROOT 

from mpl_toolkits.mplot3d import Axes3D


list_delta_L = []
list_random_point = []
list_P_vector = []

flag_long_negative = False

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

Lim_inf_theta = 0 ## grados
Lim_inf_theta_rad = np.radians(Lim_inf_theta)  ## rad

###### Dominios de variables #######
Phi = np.arange(0.001, 2 * np.pi - np.pi * 0.001, 0.001)
Radio = 10
Theta = np.arange(Lim_inf_theta_rad, (np.pi/2) - np.pi * 0.001, 0.001)    ### Semi-esfera de radio 100 unidades

##### Tamaño del plano #####
plane_side = 1.5
long_a = np.arange(-plane_side, plane_side, 0.001)
long_b = np.arange(-plane_side, plane_side, 0.001)

######### Medidas de la CCD ##########
medida_x = 1.197 / 2    # cm
medida_y = 1.587 / 2   # cm
medida_z = 0.0725   # cm

# mapeo_x = np.arange(-medida_x, medida_x, 0.01)
# mapeo_y = np.arange(-medida_y, medida_y, 0.01)
# mapeo_z = np.arange(0, medida_z, 0.0001)


mapeo_x = dimension_x(medida_x)
mapeo_y = dimension_y(medida_y)
mapeo_z = dimension_z(medida_z)

# print(mapeo_x)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.


####    Número de Puntos a Simular  ####
number_thet = 1
number_points_per_angle = 1

n_muons = number_thet * number_points_per_angle

Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

# print(mapeo_x)
#### Inicio del Bucle ####
n_muons_CCD = 0

list_delta_L_neg = []
list_delta_L_Total = []

In = datetime.datetime.now()
list_delta_L, n_muons_in_CCD, n_negative_long = func_longitud(number_thet, Theta, Theta_true, Phi, Radio, number_points_per_angle, long_a, long_b, 
                                                              medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z)
Fin = datetime.datetime.now()
print('Tiempo de cálculo para delta_L: ', Fin-In)


# if n_negative_long != 0:
#     print('Eventos con long negativa detectados: ', n_negative_long, '. Se estan volviendo a simular. ')
#     flag_long_negative = True

# if flag_long_negative:
#     while flag_long_negative:
#         list_delta_L_1, n_muons_in_CCD_1, n_negative_long = func_longitud(n_negative_long, Theta, Theta_true, Phi, Radio, 1, 
#                                                                             long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z)
    
    # if n_negative_long != 0 :
    #     print('Aun se detectaron eventos con long negativa')


#     list_delta_L_neg = list_delta_L_neg + list_delta_L_1
    
#     if n_negative_long == 0:
#         flag_long_negative = False

# list_delta_L_Total = list_delta_L + list_delta_L_neg
list_delta_L_Total = list_delta_L 
        
Final = datetime.datetime.now()

print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

print('Muones Simulados: ', n_muons)
print('Tamaño de los planos: (' + str(plane_side * 2) + ' x ' + str(plane_side * 2) + ') cm^2')
print('Muones que Impactaron la CCD: ', n_muons_in_CCD)
# print('Muones que tuvieron una longitud negativa: ', n_negative_long)

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

BINS = 450
array_Delta_L = np.array(list_delta_L_Total)
# max_DeltaL = np.max()
# print()

fig, axs = plt.subplots(figsize=[7,5])
# hist_long, bins_edges_long, _ = axs.hist(array_Delta_L, bins = BINS, label = 'Eventos Simulados: ' + str(number_thet * number_points_per_angle))
hist_long, bins_edges_long, _ = axs.hist(array_Delta_L, color = 'b', bins = BINS, label = 'Eventos Simulados: ' + str(n_muons_in_CCD),  histtype = 'step')
index_max_long =  np.argmax(hist_long)
pico =  hist_long[index_max_long]
# print(pico)
# axs.set_xlim(0, 0.5)

axs.vlines([0.0725], ymin = 0, ymax = pico, colors='k', linestyles='dashed', label = 'Grosor de la CCD: 0.0725 cm')
axs.set_xlim(0, 0.2)

axs.set_xlabel('Distancia (cm)')
axs.legend()
fig.suptitle(r'Distribución de Longitudes', y = 0.95, size = 20)
plt.show()

# file_name = 'array_Delta_L_' + str(n_muons) + '.pkl'
# file_object_histogram = open(file_name, 'wb')
# pkl.dump(array_Delta_L, file_object_histogram) ## Save the dictionary with all info 
# file_object_histogram.close()

# h = ROOT.TH1F("h", "Titutlo", 500, 0, 2)
# for element in list_delta_L:
#     h.Fill(element)


# canv = ROOT.TCanvas("canv", "Titulo 2", 700, 400)

# canv.cd()
# h.Draw()
# canv.Print("/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Simulacion_ab_initio/plot.pdf")

# print('Arrray saved in ', current_path + '/' + file_name, ' as a binary file.')
# print('To open use library "pickle". ')