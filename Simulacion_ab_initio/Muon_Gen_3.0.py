import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import *
import datetime
import os
# import pickle as pkl
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
    Radio = 4     ## cm
    # Radio = 100     ## cm


    #### Tamaño de los planos tangentes a la esfera ####
    # Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
    plane_size = 1  ## cm
    # plane_size = 75  ## cm
    long_a = np.arange(-plane_size, plane_size, 0.001)
    long_b = np.arange(-plane_size, plane_size, 0.001)


    # ######### Medidas de la CCD (para centrarla en el origen) ##########
    medida_x = 1.197 / 2    # cm
    medida_y = 1.587 / 2   # cm
    medida_z = 0.0725   # cm

    # medida_x = 100 / 2    # cm
    # medida_y = 10 / 2   # cm
    # medida_z = 10   # cm


    #### Arreglos de los valores para mapear la CCD ####
    mapeo_x = dimension_x(medida_x)
    mapeo_y = dimension_y(medida_y)
    mapeo_z = dimension_z(medida_z)

    # Max_energy = 1000000 # En MeV
    # Energy = np.arange(1, Max_energy) # En MeV

    Max_energy = 100 # En GeV (no pasar de 100 GeV o la simulación se hará mas lenta)
    Energy = np.arange(1, Max_energy, 0.001) # En GeV

    #### Distribución angular de theta (Distribución angular de Smith-Duller) ####
    Theta_true = dis_angular(Theta) 

    ### Número de muones a simular ### 
    number_thet = 100      ## Valores de un ángulo Theta.

    
    number_points_per_angle = 1  ## Valores aleatorios sobre cada plano.
    n_muons = number_thet * number_points_per_angle ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    # print(os.environ)
    
    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, nmuons_in_CCD  = muon_generator_3(Energy, number_thet=n_muons, Theta=Theta, Theta_true=Theta_true, Phi=Phi, Radio=Radio, 
                                number_points_per_angle=number_points_per_angle,long_a=long_a, long_b=long_b, medida_x=medida_x, 
                                medida_y=medida_y, medida_z=medida_z, mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z)


    #### TTree file in CCD ###
    N_Muons = array('f', [-9999])
    Thet_Rad = array('f', [-9999])
    Thet_Deg = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Phi_Deg = array('f', [-9999])
    Energy_array = array('f', [-9999])
    DeltaL_array = array('f', [-9999])
    Energy_Landau_array = array('f', [-9999])

    Thet_Deg_pri = array('f', [-9999])
    Phi_Deg_pri = array('f', [-9999])
    Energy_pri = array('f', [-9999])

    file_root_name = 'Sim_ab_initio_NMUONS_' + str(number_thet) + '.root'
    # file_root_name = 'Sim_ab_initio_Barra_NMUONS_' + str(number_thet) + '.root'
    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')

    tree.Branch('nmuon',N_Muons, 'nmuon/F' )
    tree.Branch('thet', Thet_Rad, 'thet/F')
    tree.Branch('phi', Phi_Rad, 'phi/F')
    tree.Branch('epri', Energy_array, 'epri/F')
    tree.Branch('l', DeltaL_array, 'l/F')
    tree.Branch('edep', Energy_Landau_array, 'edep/F')

    for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
        N_Muons[0] = dict_muons['NMuon'][i]
        Thet_Rad[0] = dict_muons['Theta(Rad)'][i]
        # print(Thet_Deg[0])
        #print(f'Ei={i} Energy_Landau={dict_muons_in_CCD}') 
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        DeltaL_array[0] = dict_muons['Delta_L(cm)'][i]
        Energy_Landau_array[0] = dict_muons['Energy_Landau(KeV)'][i]
        # print(Energy_Landau_array[0])
        # th_deg = dict_muons['Theta(Deg)'][0]
        tree.Fill()

    tree.Write()

    # Fin = datetime.datetime.now()
    # print('Tiempo de cálculo para hacer el tree_ROOT: ', Fin-In, end='\n\n')

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    print('Muones que impactaron en la CCD: ', nmuons_in_CCD)
    # print('TTree primary file saved in ' + file_root_name_1)
    print('TTree muons in CCd file saved in ' + file_root_name)

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


