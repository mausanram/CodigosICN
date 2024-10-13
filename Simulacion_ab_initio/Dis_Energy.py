import datetime

# In = datetime.datetime.now()
import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
# from random_geometry_points.plane import Plane   ### Para instalar utilizar "pip install random-geometry-points"
from funciones_Sim_ab_initio import *
import os
import pickle as pkl
import scipy.integrate as integrate
import pickle

# Fin = datetime.datetime.now()
# print('Tiempo de cálculo para importar librerías: ', Fin-In, end='\n\n')

### Configuración de estilo de las gráficas ###
# plt.rcParams.update({
#     "image.origin": "lower",
#     "image.aspect": 1,
#     #"text.usetex": True,
#     "grid.alpha": .5,
#     "axes.linewidth":2,
#     "lines.linewidth" : 1,
#     "font.size":    15.0,
#     "xaxis.labellocation": 'right',  # alignment of the xaxis label: {left, right, center}
#     "yaxis.labellocation": 'top',  # alignment of the yaxis label: {bottom, top, center}
#     "xtick.top":           True ,  # draw ticks on the top side
#     "xtick.major.size":    8    ,# major tick size in points
#     "xtick.minor.size":    4      ,# minor tick size in points
#     "xtick.direction":     'in',
#     "xtick.minor.visible": True,
#     "ytick.right":           True ,  # draw ticks on the top side
#     "ytick.major.size":    8    ,# major tick size in points
#     "ytick.minor.size":    4      ,# minor tick size in points
#     "ytick.direction":     'in',
#     "ytick.minor.visible": True,
#     "ytick.major.width":   2   , # major tick width in points
#     "ytick.minor.width":   1 ,
#     "xtick.major.width":   2   , # major tick width in points
#     "xtick.minor.width":   1 ,
#     "legend.framealpha": 0 ,
#     "legend.loc": 'best',
# })

current_path = os.getcwd()

def main():
    # list_rand_thet = []
    # list_rand_phi = []

    list_rand_thet_deg = []
    list_rand_phi_deg = []

    list_random_energy = []
    list_energy_Landau = []

    muon_in_bucle = 0

    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)

    ##### Cortes en el ángulo theta ####
    Lim_inf_theta_deg = 0 # grados
    Lim_inf_theta_rad = np.radians(Lim_inf_theta_deg)  ## rad

    ###### Rango de los ángulos theta y phi  #######
    Phi = np.arange(0, 2 * np.pi, 0.001) ## rad
    Theta = np.arange(Lim_inf_theta_rad, np.pi/2, 0.001) ## rad

    #### Mapeo de la energía (se hace en escala logarítmica para tener valores igualmente distribuidos) ####
    # E_in = 10 ** (-2)   ### Límite inferior
    # E_fin = 10 ** 8     ### Límite superior
    # N = 1000    ### Número de puntos

    Max_energy = 1000000 # En MeV
    Energy = np.arange(1, Max_energy) # En MeV
    # Energy = np.linspace(1, Max_energy, Max_energy)

    Max_energy = 100 # En MeV
    Energy = np.arange(1, Max_energy, 0.001) # En MeV

    #### Distribución angular de theta (Distribución angular de Smith-Duller) ####
    Theta_true = dis_angular(Theta) 

    ### Número de muones a simular ### 
    number_thet = 1000    ## Valores de un ángulo Theta.
    number_points_per_angle = 1  ## Valores aleatorios sobre cada plano.
    n_muons = number_thet * number_points_per_angle ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    for i in np.arange(0,number_thet):
        # In = datetime.datetime.now()
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Random_th: ', Fin-In)

        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes
        # print(Random_th[0])
        Random_th_deg = np.degrees(Random_th[0]) ## El ángulo theta en grados
        Random_phi_deg = np.degrees(Random_phi) ## El ángulo phi en grados

        # list_rand_thet.append(Random_th[0])
        # list_rand_phi.append(Random_phi)
        # list_rand_thet_deg.append(Random_th_deg)
        # list_rand_phi_deg.append(Random_phi_deg)

        # list_dis_Energy = []
        In = datetime.datetime.now()
        # for energy in Energy:   ## Aquí se crea la distribución de Smith-Duller en MeV
        dis_Energy = dis_energy(Energy, Random_th[0], units=1)

        # print(dis_Energy[0:20])

            # list_dis_Energy.append(dis_Energy)
        Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para distribucion_En: ', Fin-In)

        # In = datetime.datetime.now()
        # Random_energy = rand.choices(Energy, list_dis_Energy) ## Escoje una energía segun la distribución de Smith-Duller 
        Random_energy = np.around(rand.choices(Energy, dis_Energy)[0],3) * 1000 ## Escoje una energía segun la distribución de Smith-Duller (GeV to MeV)


        # print('Energy_pri: ', Random_energy, ' MeV')
        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Random_energy: ', Fin-In)

        os.environ["EN_SMITH"] = str(Random_energy)

        ## Para la laptop en el ICN  ##
        # new_env = subprocess.run(["root", "-l", "-b", "/home/labdet/Documents/MauSan/Programas/Repositorio_Git/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"],
        #                      capture_output=True)

        In = datetime.datetime.now()
        ## Para la computadora de casa ##
        new_env = subprocess.run(["root", "-l", "-b", "-n", "/home/bruce/Documents/Programas/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
                                    capture_output=True)
        Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para new_env: ', Fin-In, end='\n\n')

        ## Para el CLUSTER ##
        # new_env = subprocess.run(["root", "-l", "-b", "/home/icn/mausanram/Software/CodigosICN/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
        #                             capture_output=True)


        Random_energy_Landau = float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1])
        # print(Random_energy_Landau)


        # list_rand_thet.append(Random_th[0])
        # list_rand_phi.append(Random_phi)
        list_rand_thet_deg.append(Random_th_deg)
        list_rand_phi_deg.append(Random_phi_deg)
        list_random_energy.append(Random_energy)
        list_energy_Landau.append(Random_energy_Landau)

        muon_in_bucle += 1

        print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')
        
    # print(dis_Energy[0:20])

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    dict_to_save_pkl = {'Theta_Deg' : list_rand_thet_deg, 'Phi_Deg' : list_rand_phi_deg, 'Energy_SM' : list_random_energy, 'Energy_DP' : list_energy_Landau}

    file_name = 'Simulacion_ab_initio_MuonesSim' + str(number_thet) + '.pkl'

    file_object_histogram = open(file_name, 'wb')
    pickle.dump(dict_to_save_pkl, file_object_histogram) ## Save the dictionary with all info 
    file_object_histogram.close()

    print('Dictionary saved in', current_path + '/' + file_name, ' as a binary file. To open use library "pickle". ')


    # fig, axs_all_angle = plt.subplots(1,2, figsize = [15, 10])
    # axs_all_angle[0].hist(list_rand_thet_deg, bins = 40,  histtype = 'step')
    # axs_all_angle[0].set_xlabel(r'Angle (°)')
    # axs_all_angle[0].set_title(r'Angular $\theta$ distribution')


    # axs_all_angle[1].hist(list_rand_phi_deg, bins = 40,  histtype = 'step')
    # axs_all_angle[1].set_xlabel(r'Angle (°)')
    # axs_all_angle[1].set_title(r'Angular $\phi$ distribution')

    # plt.show()

    # fig, axs_all_energy = plt.subplots(1,2, figsize = [15, 10])
    # axs_all_energy[0].hist(list_random_energy, bins = 100,  histtype = 'step')
    # axs_all_energy[0].set_xlabel(r'Energy (MeV)')
    # axs_all_energy[0].set_title(r'Energy distribution')

    # axs_all_energy[1].hist(list_energy_Landau, bins = 100,  histtype = 'step')
    # axs_all_energy[1].set_xlabel(r'Energy (MeV)')
    # axs_all_energy[1].set_title(r'Energy DP distribution')


    # plt.show()

    
if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)
