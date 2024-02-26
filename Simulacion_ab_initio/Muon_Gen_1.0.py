import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import *
import datetime
import os
import pickle as pkl
import ROOT 

current_path = os.getcwd()
Inicio = datetime.datetime.now()
print('Hora de inicio del cálculo: ', Inicio)

###### Rango de los ángulos theta y phi  #######
Phi = np.arange(0, 2 * np.pi, 0.001) ## rad
Theta = np.arange(0, np.pi/2, 0.001) ## rad

##### Valor del radio de la semi-esfera #####
Radio = 100     ## cm


#### Tamaño de los planos tangentes a la esfera ####
# Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
plane_size = 1  ## cm
long_a = np.arange(-plane_size, plane_size, 0.001)
long_b = np.arange(-plane_size, plane_size, 0.001)


# ######### Medidas de la CCD (para centrarla en el origen) ##########
# medida_x = 1.197 / 2    # cm
# medida_y = 1.587 / 2   # cm
# medida_z = 0.0725   # cm


# #### Arreglos de los valores para mapear la CCD ####
# mapeo_x = dimension_x(medida_x)
# mapeo_y = dimension_y(medida_y)
# mapeo_z = dimension_z(medida_z)


#### Mapeo de la energía (se hace en escala logarítmica para tener valores igualmente distribuidos) ####
E_in = 10 ** (-2)   ### Límite inferior
E_fin = 10 ** 8     ### Límite superior
N = 1000    ### Número de puntos
Energy = Energy_list(E_in, E_fin, N)

#### Distribución angular de theta (Distribución angular de Smith-Duller) ####
Theta_true = dis_angular(Theta) 

### Número de muones a simular ### 
number_thet = 10    ## Valores de un ángulo Theta.
number_points_per_angle = 10  ## Valores aleatorios sobre cada plano.
n_muons = number_thet * number_points_per_angle ## Número total de muones que se simularán.

print('Se simularán ' + str(n_muons) + ' muones.')

## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
dict_muons = muon_generator_1(Radio, long_a, long_b, number_thet, number_points_per_angle, Theta, Theta_true, Phi, Energy)

# muons_dataFrame = pd.DataFrame(dict_muons)

# muons_dataFrame.to_csv('muons_data.txt', sep='\t')

tree = ROOT.TTree('T', 'Mi primer arbol')
# print(type(tree))
tree.Branch('Thet_Rad', dict_muons['Theta(Rad)'][0], 'Theta(Rad)/F')
tree.Branch('Thet_Deg', dict_muons['Theta(Deg)'][0], 'Theta(Deg)/F')
# print(tree)

for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
    th_rad = dict_muons['Theta(Rad)'][i]
    th_deg = dict_muons['Theta(Deg)'][i]

    tree.Fill()

# for element in dict_muons['Theta(Rad)']:
#     ThetRad.Fill(element)

#     # tree.Fill()

# tree.Branch('Theta(Deg)', dict_muons['Theta(Deg)'][0], 'Theta(Deg)/F')
# for element in dict_muons['Theta(Deg)']:
#     tree.Fill()

# tree.Branch('Phi(Rad)', dict_muons['Phi(Rad)'][0], 'Phi(Rad)/F')
# for element in dict_muons['Phi(Rad)']:
#     tree.Fill()

tree.Print()
# tree.Draw()
# print(tree.GetBranch('Theta(Rad)').GetEntries())
# tree.Write()

Final = datetime.datetime.now()
print('Hora final de cálculo: ', Final)
print('Tiempo de cálculo: ', Final-Inicio)

# fig, axs = plt.subplots(figsize=[7,5])
# # axs.plot(Theta, 70 * Theta_true)
# axs.hist(np.array(dict_muons['Theta(Rad)']), bins = 110)
# fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
# plt.show()

# fig, axs = plt.subplots(figsize=[7,5])
# axs.hist(np.array(dict_muons['Energy(MeV)']), bins = 110)
# fig.suptitle(r'Distribución de la Energía', y = 0.95, size = 20)
# plt.show()


