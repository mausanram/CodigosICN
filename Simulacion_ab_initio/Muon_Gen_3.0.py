print('Se están importando todas las paqueterías')
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

    pixel_size = 0.0015 # cm

    ### ======== Dimensiones originales ============ ### (CHECAR LAS DIMENSIONES DEl RADIO Y DE LOS PLANOS)
    # sizex_pixels = 1058 # px
    # sizey_pixels = 1278 # px (Esta dimensión debe ser mas larga que X)
    ### ============================================ ###

    # ### ======== Dimensiones Nuevas ================ ###
    sizex_pixels = 250 # px
    sizey_pixels = 529 # px (Esta dimensión debe ser mas larga que X)
    # ### ============================================ ###

    ### =============== Valor del radio de la semi-esfera =========== ###
    # Radio = 12     ## cm
    Radio = 8     ## cm
    ### ============================================================= ###
    
    ### =========== Tamaño de los planos tangentes a la esfera ============ ###
    # Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
    # half_plane_size = 1.2  ## cm
    half_plane_size = 0.75  ## cm
    ### ==================================================================== ###


    # ######### Medidas de la CCD (para centrarla en el origen) ##########
    medida_x = np.around((sizex_pixels * pixel_size) / 2, 4)    # cm
    medida_y = np.around((sizey_pixels * pixel_size)  / 2, 4)   # cm
    medida_z = 0.0725 / 2  # cm

    # print(medida_x, medida_y, medida_z)

    # medida_x = 1.587 / 2   # cm
    # medida_y = 1.917 / 2 # cm
    # medida_z = 0.0725 / 2  # cm

    #### Arreglos de los valores para mapear la CCD ####
    stepxy = 0.0001
    stepz = 0.00001
    mapeo_x = dimension_x(medida_x, stepxy)
    mapeo_y = dimension_y(medida_y, stepxy)
    mapeo_z = dimension_z(medida_z, stepz)

    print(mapeo_x[-1], mapeo_x[0])
    print(mapeo_y[-1], mapeo_y[0])
    print(mapeo_z[-1], mapeo_z[0])

    ### Número de muones a simular ### 
    n_muons = 1000000

    # niterations = n_muons / nmuons_perbucle

    # print('Se simularán ' + str(n_muons) + ' muones y se guardarán en objetos TTree de ' + str(nmuons_perbucle) + ' elementos.')
    print('Se simularán ' + str(n_muons))

    # for iteration in np.arange(0, niterations):
    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, nmuons_in_CCD  = muon_generator_3(number_thet=n_muons, Radio=Radio,  
                                                medida_x=medida_x, medida_y=medida_y, medida_z=medida_z, 
                                                mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z, 
                                                half_plane_size=half_plane_size)

    list_nmuons = np.arange(0, len(dict_muons['Theta(Rad)']))

    #### TTree file in CCD ###
    N_Muons = array('f', [-9999])
    Thet_Rad = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Energy_array = array('f', [-9999])
    DeltaL_array = array('f', [-9999])
    Energy_Landau_array = array('f', [-9999])
    Kappa_array = array('f', [-9999])

    # file_root_name = '/home/bruce/Documents/Programas/Simulacion_ab_initio/treesROOT_CCD/Sim_ab_initio_NMUONS_' + str(number_thet) + '_PLANES_' + str(half_plane_size * 2) +'x' + str(half_plane_size * 2) + '_RADIO_' + str(Radio) + '_' + str(iteration) + '.root'
    # file_direction = '/home/bruce/Documents/Programas/Simulacion_ab_initio/treesROOT_CCD/10k/'
    # file_root_name = 'Sim_ab_initio_NMUONS_' + str(nmuons_perbucle) + '_PLANES_' + str(half_plane_size * 2) +'x' + str(half_plane_size * 2) + '_RADIO_' + str(Radio) + '_' + str(int(iteration)) + '.root'
    # file_root_name = 'Sim_ab_initio_NMUONS_' + str(nmuons_perbucle) + '_PLANES_' + str(int(half_plane_size * 2)) +'x' + str(int(half_plane_size * 2)) + '_RADIO_' + str(Radio) + '_CCDSIZE_' + str(int(sizex_pixels))+ 'x' + str(int(sizey_pixels))+'_.root'
    file_root_name = 'Sim_ab_initio_NMUONS_' + str(n_muons) + '_PLANES_' + str(half_plane_size * 2) +'x' + str(half_plane_size * 2) + '_RADIO_' + str(Radio) + '_CCDSIZE_' + str(int(sizex_pixels))+ 'x' + str(int(sizey_pixels))+ '_SIGMA_LV_1.0' + '_.root'

    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')

    tree.Branch('nmuon',N_Muons, 'nmuon/F' )
    tree.Branch('thet', Thet_Rad, 'thet/F')
    tree.Branch('phi', Phi_Rad, 'phi/F')
    tree.Branch('epri', Energy_array, 'epri/F')
    tree.Branch('l', DeltaL_array, 'l/F')
    tree.Branch('edep', Energy_Landau_array, 'edep/F')
    tree.Branch('kappa', Kappa_array, 'kappa/F')

    print("\nLong list thet: ", len(dict_muons['Theta(Rad)']))
    print("Long list phi: ", len(dict_muons['Phi(Rad)']))
    print("Long list SD: ", len(dict_muons['Energy-SD(MeV)']))
    print("Long list L: ", len(dict_muons['Delta_L(cm)']))
    print("Long list Edep: ", len(dict_muons['Energy_Landau(KeV)']))
    print("Long list Edep: ", len(dict_muons['kappa']))

    for i in np.arange(0, len(dict_muons['Delta_L(cm)'])):
        N_Muons[0] = list_nmuons[i]
        Thet_Rad[0] = dict_muons['Theta(Rad)'][i]
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        DeltaL_array[0] = dict_muons['Delta_L(cm)'][i]
        Energy_Landau_array[0] = dict_muons['Energy_Landau(KeV)'][i]
        Kappa_array[0] = dict_muons['kappa'][i]

        # print(Thet_Rad[0], Phi_Rad[0], Energy_array[0], Energy_Landau_array[0])

        tree.Fill()

    tree.Write()
    file.Close()

    del tree, N_Muons,Thet_Rad, Phi_Rad, Energy_array, DeltaL_array, Energy_Landau_array, dict_muons, list_nmuons

    # Fin = datetime.datetime.now()
    # print('Tiempo de cálculo para hacer el tree_ROOT: ', Fin-In, end='\n\n')

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    print('Muones que impactaron en la CCD: ', nmuons_in_CCD)
    print('TTree primary file saved in ' + file_root_name)
    # print('All TTree muons in CCd file saved in /home/bruce/Documents/Programas/Simulacion_ab_initio/treesROOT_CCD')

if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)


