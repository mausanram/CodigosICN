import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import *
import datetime
import os
import pickle as pkl
# import ROOT 
from ROOT import TFile, TTree
from array import array


def main():
    current_path = os.getcwd()
    Inicio = datetime.datetime.now()
    print('Hora de inicio del cálculo: ', Inicio)


    ##### Cortes en el ángulo theta ####
    Lim_inf_theta_deg = 0 # grados
    Lim_inf_theta_rad = np.radians(Lim_inf_theta_deg)  ## rad

    ###### Rango de los ángulos theta y phi  #######
    Phi = np.arange(0, 2 * np.pi, 0.001) ## rad
    Theta = np.arange(Lim_inf_theta_rad, np.pi/2, 0.001) ## rad

    ##### Valor del radio de la semi-esfera #####
    Radio = 100     ## cm


    #### Tamaño de los planos tangentes a la esfera ####
    # Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
    plane_size = 1.5  ## cm
    long_a = np.arange(-plane_size, plane_size, 0.001)
    long_b = np.arange(-plane_size, plane_size, 0.001)


    # ######### Medidas de la CCD (para centrarla en el origen) ##########
    medida_x = 1.197 / 2    # cm
    medida_y = 1.587 / 2   # cm
    medida_z = 0.0725   # cm


    #### Arreglos de los valores para mapear la CCD ####
    mapeo_x = dimension_x(medida_x)
    mapeo_y = dimension_y(medida_y)
    mapeo_z = dimension_z(medida_z)


    #### Mapeo de la energía (se hace en escala logarítmica para tener valores igualmente distribuidos) ####
    # E_in = 10 ** (-2)   ### Límite inferior
    # E_fin = 10 ** 8     ### Límite superior
    # N = 1000    ### Número de puntos
    # Energy = Energy_list(E_in, E_fin, N)

    # Max_energy = 1000000 # E MeV
    # Energy = np.arange(1, Max_energy) # En MeV

    Max_energy = 100 # En GeV
    Energy = np.arange(1, Max_energy, 0.001) # En GeV

    #### Distribución angular de theta (Distribución angular de Smith-Duller) ####
    Theta_true = dis_angular(Theta) 

    ### Número de muones a simular ### 
    number_thet = 1    ## Valores de un ángulo Theta.

    
    number_points_per_angle = 1  ## Valores aleatorios sobre cada plano.
    n_muons = number_thet * number_points_per_angle ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    # print(os.environ)
    
    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, _, _  = muon_generator_3(Energy, number_thet=n_muons, Theta=Theta, Theta_true=Theta_true, Phi=Phi, Radio=Radio, 
                                number_points_per_angle=number_points_per_angle,long_a=long_a, long_b=long_b, medida_x=medida_x, 
                                medida_y=medida_y, medida_z=medida_z, mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z)



    # print(os.environ)

    # muons_dataFrame = pd.DataFrame(dict_muons)
    # muons_dataFrame.to_csv('muons_data.txt', sep='\t')
    In = datetime.datetime.now()

    file_root_name = 'Sim_ab_initio_NMUONS_' + str(number_thet) + '.root'
    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')
    # print(type(tree))

    # print('Longitud la lista de "Energy_Landau": ', len(dict_muons['Energy_Landau']))
    # print('Longitud la lista de "Theta(Rad)": ', len(dict_muons['Theta(Rad)']))

    # print(dict_muons['Energy_Landau'])
    # print(dict_muons['Theta(Deg)'])

    Thet_Rad = array('f', [-9999])
    Thet_Deg = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Phi_Deg = array('f', [-9999])
    Energy_array = array('f', [-9999])
    DeltaL_array = array('f', [-9999])
    Energy_Landau_array = array('f', [-9999])

    tree.Branch('Thet_Rad', Thet_Rad, 'Thet_Rad/F')
    tree.Branch('Thet_Deg', Thet_Deg, 'Thet_Deg/F')
    tree.Branch('Phi_Rad', Phi_Rad, 'Phi_Rad/F')
    tree.Branch('Phi_Deg', Phi_Deg, 'Phi_Deg/F')
    tree.Branch('Energy', Energy_array, 'Energy/F')
    tree.Branch('DeltaL', DeltaL_array, 'DeltaL/F')
    tree.Branch('Energy_DP', Energy_Landau_array, 'Energy_DP/F')


    for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
        Thet_Rad = dict_muons['Theta(Rad)'][i]
        Thet_Deg[0] = dict_muons['Theta(Deg)'][i]
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Phi_Deg[0] = dict_muons['Phi(Deg)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        DeltaL_array[0] = dict_muons['Delta_L'][i]
        Energy_Landau_array[0] = dict_muons['Energy_Landau'][i]
        # th_deg = dict_muons['Theta(Deg)'][0]
        tree.Fill()

    # tree.Show(-1)
    # tree.Print()
    # tree.Draw()
    # print(tree.GetBranch('Theta(Rad)').GetEntries())
    tree.Write()
    Fin = datetime.datetime.now()
    print('Tiempo de cálculo para hacer el tree_ROOT: ', Fin-In, end='\n\n')

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    print('Muones que impactaron en la CCD: ', len(dict_muons['Energy_Landau']))
    print('TTree file saved in ' + file_root_name)

    # fig, axs = plt.subplots(figsize=[7,5])
    # # axs.plot(Theta, 70 * Theta_true)
    # axs.hist(np.array(dict_muons['Theta(Rad)']), bins = 110)
    # fig.suptitle(r'Distribución angular $\theta$', y = 0.95, size = 20)
    # plt.show()

    # fig, axs = plt.subplots(figsize=[7,5])
    # axs.hist(np.array(dict_muons['Energy(MeV)']), bins = 110)
    # fig.suptitle(r'Distribución de la Energía', y = 0.95, size = 20)
    # plt.show()

if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)


