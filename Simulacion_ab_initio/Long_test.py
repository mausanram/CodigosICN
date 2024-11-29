import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
from funciones_Sim_ab_initio import *

def main():
    current_path = os.getcwd()
    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)

    ##### Valor del radio de la semi-esfera #####
    Radio = 12     ## cm

    #### Tamaño de los planos tangentes a la esfera ####
    # Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
    half_plane_size = 1.5  ## cm

    # ######### Medidas de la CCD (para centrarla en el origen) ##########
    medida_x = 1.197 / 2    # cm
    medida_y = 1.587 / 2   # cm
    medida_z = 0.0725   # cm

    #### Arreglos de los valores para mapear la CCD ####
    mapeo_x = dimension_x(medida_x)
    mapeo_y = dimension_y(medida_y)
    mapeo_z = dimension_z(medida_z)

    ### Número de muones a simular ### 
    number_thet = 1000
    n_muons = number_thet ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, nmuons_in_CCD  = muon_generator_2(number_thet=n_muons, Radio=Radio,
                                       medida_x=medida_x, medida_y=medida_y, medida_z=medida_z, 
                                       mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z, 
                                       half_plane_size = half_plane_size)


    print('Muones que impactaron en la CCD: ', nmuons_in_CCD)
    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)
    # plt.show()

if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)

# ### =========== Cara inferior =============== ###
# X1 = [-medida_x, medida_x, medida_x, -medida_x, -medida_x, medida_x]
# Y1 = [-medida_y, -medida_y, medida_y, medida_y, -medida_y, medida_y]
# Z1 = [0, 0, 0, 0, 0, 0]
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(xs=X1, ys=Y1, zs=Z1, color='k')

# ### =========== Cara superior =============== ###
# X1 = [-medida_x, medida_x, medida_x, -medida_x, -medida_x]
# Y1 = [-medida_y, -medida_y, medida_y, medida_y, -medida_y]
# Z1 = [medida_z, medida_z, medida_z, medida_z, medida_z]
# ax.plot(xs=X1, ys=Y1, zs=Z1, color='k')

# ### ========== Uniones ==================== ###
# x1= [-medida_x, -medida_x]
# y1 = [-medida_y, -medida_y]
# z1 = [0, medida_z]
# ax.plot(xs=x1, ys=y1, zs=z1, color='k')

# x1= [medida_x, medida_x]
# y1 = [-medida_y, -medida_y]
# z1 = [0, medida_z]
# ax.plot(xs=x1, ys=y1, zs=z1, color='k')

# x1= [medida_x, medida_x]
# y1 = [medida_y, medida_y]
# z1 = [0, medida_z]
# ax.plot(xs=x1, ys=y1, zs=z1, color='k')

# x1= [-medida_x, -medida_x]
# y1 = [medida_y, medida_y]
# z1 = [0, medida_z]
# ax.plot(xs=x1, ys=y1, zs=z1, color='k')


# ### ============= Proyecion del domo ================ ###
# # theta = np.linspace(0, 4 * np.pi, 100)
# # # z = np.linspace(-2, 2, 100)
# # z = 0
# # r = 12
# # x = r * np.sin(theta)
# # y = r * np.cos(theta)
# # ax.plot(x, y, z, label='parametric curve')

# # ax.scatter([0,0], [0,0], [0,6])

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# plt.show()
# ax.set_xscale()
# ax.scatter(xs=1,ys=1,zs=1)