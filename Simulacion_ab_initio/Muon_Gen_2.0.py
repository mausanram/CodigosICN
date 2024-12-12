import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import *
import datetime
import os
# import pickle as pkl

from ROOT import TFile, TTree
from array import array


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
    number_thet = 500000    
    n_muons = number_thet ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    # print(os.environ)
    
    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, nmuons_in_CCD  = muon_generator_2(number_thet=n_muons, Radio=Radio,
                                       medida_x=medida_x, medida_y=medida_y, medida_z=medida_z, 
                                       mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z, 
                                       half_plane_size = half_plane_size)

    # print('Numero de muones que impactaron: ', nmuons_in_CCD)
    list_nmuons = np.arange(0, len(dict_muons['Theta(Rad)']))

    #### TTree file in CCD ###
    N_Muons = array('f', [-9999])
    Thet_Rad = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Energy_array = array('f', [-9999])
    DeltaL_array = array('f', [-9999])
    Energy_Landau_array = array('f', [-9999])

    file_root_name = 'Sim_ab_initio_NMUONS_' + str(number_thet) + '_PLANES_' + str(half_plane_size * 2) +'x' + str(half_plane_size * 2) + '_RADIO_' + str(Radio) + '_0.root'
    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')

    tree.Branch('nmuon',N_Muons, 'nmuon/F' )
    tree.Branch('thet', Thet_Rad, 'thet/F')
    tree.Branch('phi', Phi_Rad, 'phi/F')
    tree.Branch('epri', Energy_array, 'epri/F')
    tree.Branch('l', DeltaL_array, 'l/F')
    tree.Branch('edep', Energy_Landau_array, 'edep/F')

    for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
        N_Muons[0] = list_nmuons[i]
        # N_Muons[0] = dict_muons['NMuon'][i]
        Thet_Rad[0] = dict_muons['Theta(Rad)'][i]
        # print(Thet_Deg[0])
        #print(f'Ei={i} Energy_Landau={dict_muons_in_CCD}') 
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        DeltaL_array[0] = dict_muons['Delta_L(cm)'][i]
        # Energy_Landau_array[0] = dict_muons['Energy_Landau(KeV)'][i]
        # print(Energy_Landau_array[0])
        # th_deg = dict_muons['Theta(Deg)'][0]
        tree.Fill()

    tree.Write()

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    print('Muones que impactaron en la CCD: ', nmuons_in_CCD)
    # print('TTree primary file saved in ' + file_root_name_1)
    print('TTree muons in CCd file saved in ' + file_root_name)

if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)


