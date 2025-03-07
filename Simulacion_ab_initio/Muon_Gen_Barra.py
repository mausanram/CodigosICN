print('Se están importando todas las paqueterías')
import numpy as np
import mpmath as mp
import random as rand
import pandas as pd
import matplotlib.pyplot as plt
from funciones_Sim_ab_initio import dimension_x_barr, dimension_y_barr, dimension_z_barr, muon_generator_BARRA
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

    ##### Valor del radio de la semi-esfera #####
    # Radio = 8     ## cm
    Radio = 450     ## cm

    #### Tamaño de los planos tangentes a la esfera ####
    # Son planos simétricos de tamaño (2 * plane_side x 2 * plane_side)
    # plane_size = 2  ## cm
    half_plane_size = 75  ## cm

    # ######### Medidas de la CCD (para centrarla en el origen) ##########
    # medida_x = 1.197 / 2    # cm
    # medida_y = 1.587 / 2   # cm
    # medida_z = 0.0725   # cm

    medida_x = 10 / 2    # cm
    medida_y = 100 / 2   # cm
    medida_z = 10 / 2  # cm


    #### Arreglos de los valores para mapear la CCD ####
    mapeo_x = dimension_x_barr(medida_x)
    mapeo_y = dimension_y_barr(medida_y)
    mapeo_z = dimension_z_barr(medida_z)

    ### Número de muones a simular ### 
    number_thet = 100000      ## Valores de un ángulo Theta.

    n_muons = number_thet  ## Número total de muones que se simularán.

    print('Se simularán ' + str(n_muons) + ' muones.')

    # print(os.environ)
    
    ## Se simulan los muones, se genera un diccionario con la información de cada evento (Theta, Phi, Energía) ##
    dict_muons, nmuons_in_CCD  = muon_generator_BARRA(number_thet=n_muons, Radio=Radio, 
                                                      medida_x=medida_x, medida_y=medida_y, medida_z=medida_z, 
                                                      mapeo_x=mapeo_x, mapeo_y=mapeo_y, mapeo_z=mapeo_z, 
                                                      half_plane_size =half_plane_size)

    list_nmuons = np.arange(0, len(dict_muons['Theta(Rad)']))

    #### TTree file in CCD ###
    N_Muons = array('f', [-9999])
    Thet_Rad = array('f', [-9999])
    Phi_Rad = array('f', [-9999])
    Energy_array = array('f', [-9999])
    DeltaL_array = array('f', [-9999])
    Energy_Landau_array = array('f', [-9999])
    Kappa_array = array('f', [-9999])


    file_root_name = 'Sim_ab_initio_Barra_NMUONS_' + str(number_thet)  + '_PLANES_' + str(half_plane_size * 2) +'x' + str(half_plane_size * 2) + '_RADIO_' + str(Radio)  + '.root'
    file = TFile.Open(file_root_name, "RECREATE")
    tree = TTree('tree', 'tree')

    tree.Branch('nmuon',N_Muons, 'nmuon/F' )
    tree.Branch('thet', Thet_Rad, 'thet/F')
    tree.Branch('phi', Phi_Rad, 'phi/F')
    tree.Branch('epri', Energy_array, 'epri/F')
    tree.Branch('l', DeltaL_array, 'l/F')
    tree.Branch('edep', Energy_Landau_array, 'edep/F')
    tree.Branch('kappa', Kappa_array, 'kappa/F')


    for i in np.arange(0, len(dict_muons['Theta(Rad)'])):
        N_Muons[0] = list_nmuons[i]
        Thet_Rad[0] = dict_muons['Theta(Rad)'][i]
        Phi_Rad[0] = dict_muons['Phi(Rad)'][i]
        Energy_array[0] =  dict_muons['Energy-SD(MeV)'][i] 
        DeltaL_array[0] = dict_muons['Delta_L(cm)'][i]
        Energy_Landau_array[0] = dict_muons['Energy_Landau(KeV)'][i]
        Kappa_array[0] = dict_muons['kappa'][i]


        tree.Fill()

    tree.Write()
    file.Close()

    # Fin = datetime.datetime.now()
    # print('Tiempo de cálculo para hacer el tree_ROOT: ', Fin-In, end='\n\n')

    Final = datetime.datetime.now()
    print('Hora final de cálculo: ', Final)
    print('Tiempo de cálculo: ', Final-Inicio)

    print('Muones que impactaron en la CCD: ', nmuons_in_CCD)
    # print('TTree primary file saved in ' + file_root_name_1)
    print('TTree muons in CCd file saved in ' + file_root_name)


if __name__ == "__main__":
    exitcode = main()
    exit(code = exitcode)


