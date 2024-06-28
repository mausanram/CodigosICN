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

########## Checar que la constante de ajuste sea la correcta y se ajuste bien a la dstribución ########

Lim_inf_theta = 22 ## grados
Lim_inf_theta_rad = np.radians(Lim_inf_theta)  ## rad

###### Definiciones y dominios #######
Phi = np.arange(0, 2 * np.pi, 0.001)
Radio = 100
Theta = np.arange(Lim_inf_theta_rad, np.pi/2, 0.001)    ### Semi-esfera de radio 100 unidades

Theta_deg = []
for angle in Theta:
    Thet_deg = np.degrees(angle)
    Theta_deg.append(Thet_deg)


long_a = np.arange(-5, 5, 0.01)
long_b = np.arange(-5, 5, 0.01)

# Energy = np.arange(10, 100000, 1)

Theta_true = dis_angular(Theta) ## Distribución angular theta real.

# print('Clase de tipo: ', type(Theta_true))
# print('Longitud Total: ', len(Theta_true))
# print('Longitud: ', Theta[-1])

number_thet = 61120
number_points_per_angle = 1

n_muons = number_thet * number_points_per_angle

Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

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

Bins = 100
delta_B = Theta[-1]/ Bins
Const_Normal =  3 * number_thet * delta_B
print('La constante de ajuste de la curva es: ', np.round(Const_Normal, 3))

list_random_thet_deg = []
for angle in list_random_thet:
    angle_deg = np.degrees(angle)
    # print(angle, angle_deg)
    list_random_thet_deg.append(angle_deg)

list_random_phi_deg = []
for angle in list_random_phi:
     angle_deg = np.degrees(angle)
     list_random_phi_deg.append(angle_deg)


Final = datetime.datetime.now()
print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

# ----------------------------------------------------------------------------------------- #

### En radianes ###
########### Gráfica de Theta ###############
# fig, axs = plt.subplots(figsize=[10,5])
# axs.plot(Theta, Const_Normal * Theta_true, label = str(np.round(Const_Normal, 2)) + r'*$sin \theta cos^2 \theta$')
# axs.hist(np.array(list_random_thet), bins = Bins, range=[0, Theta[-1]],  label='Eventos Simulados: \n ' + str(n_muons))
# axs.legend()
# fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
# axs.set_xlabel('Ángulo (Rad)')
# plt.show()

########### Gráfica de Phi #################
# fig, axs = plt.subplots(figsize=[7,5])
# axs.hist(list_random_phi, bins = Bins)
# fig.suptitle(r'Distribución angular $\phi$', y = 0.95, size = 20)
# plt.show()

# ----------------------------------------------------------------------------------------- #

### En grados ### 
########### Gráfica de Theta ###############
fig, axs = plt.subplots(figsize=[10,5])
axs.plot(Theta_deg, Const_Normal * Theta_true, label = str(np.round(Const_Normal, 2)) + r'*$sin \theta cos^2 \theta$')
axs.hist(np.array(list_random_thet_deg), bins = Bins, color = 'k', label='Eventos Simulados: \n ' + str(n_muons), histtype = 'step')
axs.legend()
axs.set_xlabel('Ángulo (°)')
axs.set_xlim(0,)
fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
plt.show()

########### Gráfica de Phi #################
fig, axs = plt.subplots(figsize=[8,5])
axs.hist(list_random_phi_deg, bins = Bins)
fig.suptitle(r'Distribución angular $\phi$', y = 0.95, size = 20)
axs.set_xlabel('Ángulo (°)')
plt.show()

# ----------------------------------------------------------------------------------------- #

